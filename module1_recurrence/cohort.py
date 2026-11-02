"""Oncotype DX 대상 코호트 정의 (NCCN: TAILORx, RxPONDER 근거).

TCGA · BCNB · METABRIC · SCAN-B를 같은 표준 열로 바꾼 뒤 assign_group 하나로 나눈다.
표준 열: patient_id, hr, her2, node_group(N0 / N1 / N2+ / None), tumor_gt_05cm, postmenopausal

실행(BCNB 확인용): python -m module1_recurrence.cohort
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BCNB_CLINICAL = ROOT / "BCNB" / "BCNB" / "patient-clinical-data.xlsx"

# 검사 권고: N0 & 종양 > 0.5cm, N1(1–3개) & 폐경 후
RECOMMENDED = {"N0", "N1_post"}


def nodes_to_group(n):
    if pd.isna(n):
        return None
    return "N0" if n == 0 else "N1" if n <= 3 else "N2+"


def pn_to_group(pn):
    """AJCC pN → 림프절 구간. NCCN처럼 N1mi(미세전이)는 N0로 본다."""
    pn = str(pn)
    if pn.startswith("N0") or pn == "N1mi":
        return "N0"
    if pn.startswith("N1"):
        return "N1"
    if pn.startswith(("N2", "N3")):
        return "N2+"
    return None


def assign_group(df):
    """표준 열을 가진 DataFrame → 환자별 코호트 이름."""
    def one(r):
        if not (r.hr and r.her2 == "Negative"):
            return "not_HR+HER2-"
        if r.node_group == "N2+":
            return "N2+"
        if not r.tumor_gt_05cm:
            return "T<=0.5cm"
        if r.node_group == "N0":
            return "N0"
        if r.node_group == "N1":
            return "N1_post" if r.postmenopausal else "N1_pre"
        return "node_unknown"
    return df.apply(one, axis=1)


def tcga_her2(r):
    """IHC 양성 또는 FISH 양성 → 양성, IHC 음성 또는 (equivocal & FISH 음성) → 음성."""
    ihc, fish = r.her2_status_by_ihc, r.her2_fish_status
    if ihc == "Positive" or fish == "Positive":
        return "Positive"
    if ihc == "Negative" or (ihc == "Equivocal" and fish == "Negative"):
        return "Negative"
    return "Unknown"


def standardize_tcga(clin):
    """nationwidechildrens.org_clinical_patient_brca.txt → 표준 열.

    림프절 수(H&E 양성 개수)가 없으면 pN으로 구간을 채운다.
    pN1mi는 개수에 미세전이가 섞여 있을 수 있어서 개수 대신 pN(N0)을 쓴다.
    """
    nodes = pd.to_numeric(clin.lymph_nodes_examined_he_count, errors="coerce")
    is_n1mi = clin.ajcc_nodes_pathologic_pn == "N1mi"
    from_count = nodes.map(nodes_to_group).where(~is_n1mi)
    from_pn = clin.ajcc_nodes_pathologic_pn.map(pn_to_group)
    out = pd.DataFrame({
        "patient_id": clin.bcr_patient_barcode,
        "hr": (clin.er_status_by_ihc == "Positive") | (clin.pr_status_by_ihc == "Positive"),
        "her2": clin.apply(tcga_her2, axis=1),
        "nodes": nodes,
        "node_group": from_count.fillna(from_pn),
        "node_source": from_count.notna().map({True: "count", False: "pN"}).where(
            from_count.notna() | from_pn.notna()),
        # 크기 열이 없어서 pT로 판단: T1mi, T1a가 0.5cm 이하
        "tumor_gt_05cm": ~clin.ajcc_tumor_pathologic_pt.isin(["T1mi", "T1a"]),
        "postmenopausal": clin.menopause_status.astype(str).str.startswith("Post"),
    })
    out["group"] = assign_group(out)
    return out


def standardize_bcnb(df):
    """BCNB patient-clinical-data.xlsx → 표준 열. 폐경 정보가 없어서 50세 이상을 폐경 후로 본다."""
    nodes = df["Number of lymph node metastases"]
    out = pd.DataFrame({
        "patient_id": df["Patient ID"],
        "hr": (df.ER == "Positive") | (df.PR == "Positive"),
        "her2": df.HER2,
        "nodes": nodes,
        "node_group": nodes.map(nodes_to_group),
        "node_source": "count",
        "tumor_gt_05cm": df["Tumour Size(cm)"] > 0.5,
        "postmenopausal": df["Age(years)"] >= 50,
    })
    out["group"] = assign_group(out)
    return out


def standardize_metabric(df):
    """data/metabric/hr_her2neg.csv(cBioPortal 임상) → 표준 열.

    림프절은 양성 개수만 있다(pN1mi 구분 없음). 종양 크기는 mm, 결측이면 TCGA처럼 0.5cm 초과로 본다.
    """
    nodes = df.LYMPH_NODES_EXAMINED_POSITIVE
    out = pd.DataFrame({
        "patient_id": df.patientId,
        "hr": (df.ER_STATUS == "Positive") | (df.PR_STATUS == "Positive"),
        "her2": df.HER2_STATUS,
        "nodes": nodes,
        "node_group": nodes.map(nodes_to_group),
        "node_source": "count",
        "tumor_gt_05cm": ~(df.TUMOR_SIZE <= 5),
        "postmenopausal": df.INFERRED_MENOPAUSAL_STATE == "Post",
    })
    out["group"] = assign_group(out)
    return out


SCANB_LN = {"N0": "N0", "SubMicroMet": "N0", "1to3": "N1", "4toX": "N2+"}
SCANB_AGE = "Age (5-year range, e.g., 35(31-35), 40(36-40), 45(41-45) etc.)"


def standardize_scanb(df):
    """SCAN-B 임상 표(SCANB.9206 시트) → 표준 열.

    나이가 5세 구간 상한(50 = 46–50세)으로만 있고 폐경 정보가 없어서, 51세 이상(구간 상한 55 이상)을 폐경 후로 본다.
    """
    out = pd.DataFrame({
        "patient_id": df.Patient,
        "hr": (df.ER == "Positive") | (df.PR == "Positive"),
        "her2": df.HER2,
        "nodes": df["LN.spec"],
        "node_group": df["LN.spec"].map(SCANB_LN),
        "node_source": "LN.spec",
        "tumor_gt_05cm": ~df.pT.isin(["T1mi", "T1a", "Tis"]),
        "postmenopausal": df[SCANB_AGE] > 50,
    })
    out["group"] = assign_group(out)
    return out


if __name__ == "__main__":
    bcnb = standardize_bcnb(pd.read_excel(BCNB_CLINICAL))
    print(bcnb.group.value_counts().to_dict())
    print("HR+/HER2-:", (bcnb.group != "not_HR+HER2-").sum(),
          "| recommended:", bcnb.group.isin(RECOMMENDED).sum())
