"""TCGA 임상 · 진단 슬라이드 · RNA-seq · TCGA-CDR을 환자 바코드로 연결해 manifest.csv를 만든다.

먼저 gdc_fetch를 실행해 data/tcga/raw/를 채운다.
실행: python -m module1_recurrence.data.scripts.manifest
결과: data/tcga/manifest.csv (환자 1행)
"""
import pandas as pd

from module1_recurrence.scripts.cohort import RECOMMENDED, standardize_tcga
from module1_recurrence.data.scripts.gdc_fetch import CLINICAL_FILE_NAME, RAW

OUT = RAW.parent / "manifest.csv"
OUT_OTHER = RAW.parent / "other" / "manifest.csv"
OUT_DISTANT = RAW.parent / "distant.csv"  # TCGA-BRCA 전체 환자의 원격 전이 사건 · 일수 (dm_event, dm_days)


def load_slides():
    """환자당 진단 슬라이드(DX1, DX2 …)를 한 행으로 묶는다."""
    dx = pd.read_csv(RAW / "dx_slides.tsv", sep="\t").sort_values("file_name")
    return dx.groupby("cases.0.submitter_id").agg(
        slide_file_ids=("file_id", ";".join),
        slide_file_names=("file_name", ";".join),
        n_slides=("file_id", "size"),
        slide_gb=("file_size", lambda s: round(s.sum() / 1e9, 3)),
    )


def load_rna():
    """원발 종양 RNA-seq. 파일이 2개 이상이면 샘플 바코드 순으로 첫 번째(01A)를 쓴다."""
    rna = pd.read_csv(RAW / "rnaseq.tsv", sep="\t")
    rna = rna.sort_values("cases.0.samples.0.submitter_id")
    rna = rna.drop_duplicates("cases.0.submitter_id")
    return rna.set_index("cases.0.submitter_id")[["file_id", "cases.0.samples.0.submitter_id"]].rename(
        columns={"file_id": "rna_file_id", "cases.0.samples.0.submitter_id": "rna_sample_id"})


def load_cdr():
    cdr = pd.read_excel(RAW / "TCGA-CDR-SupplementalTableS1.xlsx", sheet_name="TCGA-CDR")
    cdr = cdr[cdr.type == "BRCA"].set_index("bcr_patient_barcode")
    return cdr[["DFI", "DFI.time", "PFI", "PFI.time"]].rename(
        columns={"DFI.time": "DFI_time", "PFI.time": "PFI_time"})


NTE_FILES = ["nationwidechildrens.org_clinical_nte_brca.txt", "nationwidechildrens.org_clinical_follow_up_v4.0_nte_brca.txt",
             "nationwidechildrens.org_clinical_follow_up_v2.1_brca.txt"]


def load_distant():
    """원격 전이(new_tumor_event_type = Distant Metastasis)까지 일수. 기록이 여러 개면 가장 이른 것.

    TCGA-CDR의 PFI는 국소 재발 · 원격 전이 · 새 원발암 · 진행 · 암 사망을 합친 것이라 원격 재발 검증에는 이 값을 따로 쓴다.
    원격 전이가 없으면 마지막 접촉(OS.time, 사망 포함)에서 중도 절단. 파일은 GDC TCGA-BRCA 임상 BCR Biotab (gdc_fetch 참고).
    """
    rows = []
    for f in NTE_FILES:
        if not (RAW / f).exists():
            continue
        d = pd.read_csv(RAW / f, sep="\t", skiprows=[1, 2])
        x = d[d.new_tumor_event_type.astype(str).str.contains("Distant")]
        rows += list(zip(x.bcr_patient_barcode, pd.to_numeric(x.new_tumor_event_dx_days_to, errors="coerce")))
    days = pd.DataFrame(rows, columns=["patient_id", "days"]).groupby("patient_id").days.min()
    cdr = pd.read_excel(RAW / "TCGA-CDR-SupplementalTableS1.xlsx", sheet_name="TCGA-CDR")
    cdr = cdr[cdr.type == "BRCA"].set_index("bcr_patient_barcode")
    out = pd.DataFrame({"last_contact_days": cdr["OS.time"]})
    out["dm_event"] = out.index.isin(days.index).astype(int)
    out["dm_days"] = days.reindex(out.index).fillna(out.last_contact_days)
    return out[["dm_event", "dm_days"]]


def build():
    # 1–2행은 긴 열 이름과 CDE_ID라서 건너뛴다
    clin = pd.read_csv(RAW / CLINICAL_FILE_NAME, sep="\t", skiprows=[1, 2])
    std = standardize_tcga(clin)
    std["tss"] = std.patient_id.str.split("-").str[1]

    eligible = (
        (clin.gender == "FEMALE")
        & ~clin.ajcc_metastasis_pathologic_pm.astype(str).str.startswith("M1")
        & (clin.history_neoadjuvant_treatment != "Yes")
    )
    m = (std[eligible.values]
         .join(load_slides(), on="patient_id")
         .join(load_rna(), on="patient_id")
         .join(load_cdr(), on="patient_id"))

    print("TCGA-BRCA", len(std), "| female, M0/MX, no neoadjuvant", len(m))
    m = m[m.group != "not_HR+HER2-"]
    print("HR+/HER2-", len(m), end="")
    m = m[m.slide_file_ids.notna()]
    print(" | + DX slide", len(m), end="")
    m = m[m.rna_file_id.notna()]
    print(" | + RNA-seq", len(m))

    m["recommended"] = m.group.isin(RECOMMENDED)
    m.to_csv(OUT, index=False)
    # 원격 재발 결과는 manifest를 바꾸지 않고 따로 저장 (H&E → 재발 평가용)
    load_distant().rename_axis("patient_id").to_csv(OUT_DISTANT)

    # H&E → 재발 점수 학습을 늘리기 위한 나머지 아형 (C2 방법 4): 여성 + DX 슬라이드 + RNA, 위 코호트 밖
    # 원격 전이 검증용: 진단 시 원격 전이(M1)인 환자는 재발이 아니라서 m1 표시를 남긴다
    other = (std[(clin.gender == "FEMALE").values].join(load_slides(), on="patient_id")
             .join(load_rna(), on="patient_id").join(load_cdr(), on="patient_id").join(load_distant(), on="patient_id"))
    other["m1"] = clin.set_index("bcr_patient_barcode").ajcc_metastasis_pathologic_pm.reindex(other.patient_id).astype(str).str.startswith("M1").to_numpy()
    other = other[other.slide_file_ids.notna() & other.rna_file_id.notna() & ~other.patient_id.isin(m.patient_id)]
    other.to_csv(OUT_OTHER, index=False)
    print(f"other subtypes (female + DX + RNA, not in cohort) {len(other)} | slides GB {other.slide_gb.sum():.0f} → {OUT_OTHER}")

    summary = m.groupby("group").agg(n=("patient_id", "size"), pfi_events=("PFI", "sum"),
                                     dfi_n=("DFI", "count"), dfi_events=("DFI", "sum"))
    print(summary.astype(int))
    rec = m[m.recommended]
    print(f"recommended {len(rec)} | PFI events {int(rec.PFI.sum())} | DFI events {int(rec.DFI.sum())}")
    print("node_source", m.node_source.value_counts(dropna=False).to_dict())
    print(f"slides {int(m.n_slides.sum())} files, {m.slide_gb.sum():.0f} GB | TSS {m.tss.nunique()}")
    print("saved", OUT)


if __name__ == "__main__":
    build()
