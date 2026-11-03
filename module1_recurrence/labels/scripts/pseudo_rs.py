"""RNA-seq에 Oncotype 21-gene 공식(Paik 2004)을 적용해 유사 RS(0–100)를 계산한다.

log2 TPM → 기준 유전자 5개 평균으로 정규화 → 0–15 척도 → 그룹 점수 → RS.
계수와 그룹 하한값은 genefu oncotypedx.R 소스와 대조했다.
RNA-seq는 RT-qPCR과 척도가 달라 RS 26 기준이 맞지 않는다 (그대로 쓰면 73%가 26 이상).
그래서 위험 구간은 점수 순위로 나눈다: TAILORx 분포처럼 하위 16% low, 상위 17% high.
실행: python -m module1_recurrence.labels.scripts.pseudo_rs
결과: data/tcga/pseudo_rs.csv (patient_id, 그룹 점수, pseudo_rs, rs_group, high_risk)
"""
import numpy as np
import pandas as pd

from module1_recurrence.data.scripts.gdc_fetch import RAW

MANIFEST = RAW.parent / "manifest.csv"
RNA_DIR = RAW.parent / "rnaseq"
OUT = RAW.parent / "pseudo_rs.csv"

PROLIFERATION = ["MKI67", "AURKA", "BIRC5", "CCNB1", "MYBL2"]
CANCER_GENES = PROLIFERATION + ["ERBB2", "GRB7", "ESR1", "PGR", "BCL2", "SCUBE2",
                                "MMP11", "CTSV", "GSTM1", "CD68", "BAG1"]
REFERENCE = ["ACTB", "GAPDH", "RPLP0", "GUSB", "TFRC"]

# Paik 2004 논문 이름 → 현재 HGNC 기호 (TCGA gene_name은 HGNC 기호)
ALIASES = {"Ki67": "MKI67", "STK15": "AURKA", "Survivin": "BIRC5", "HER2": "ERBB2",
           "ER": "ESR1", "PR": "PGR", "STMY3": "MMP11", "CTSL2": "CTSV", "GUS": "GUSB"}

LOW_PCT, HIGH_PCT = 0.16, 0.17  # TAILORx: RS 0–10 약 16%, 26 이상 약 17%


def to_hgnc(expr):
    return expr.rename(columns=ALIASES)


def load_expression(patient_ids):
    """환자 × 21개 유전자 log2(TPM + 1)."""
    rows = {}
    for pid in patient_ids:
        d = pd.read_csv(RNA_DIR / f"{pid}.tsv", sep="\t", comment="#",
                        usecols=["gene_name", "tpm_unstranded"])
        rows[pid] = d[d.gene_name.isin(CANCER_GENES + REFERENCE)].groupby("gene_name").tpm_unstranded.sum()
    return np.log2(pd.DataFrame(rows).T + 1)


def normalize(log_expr):
    """기준 유전자 평균을 빼고, 유전자마다 코호트 안에서 0–15로 min-max (genefu do.scaling과 같음)."""
    x = log_expr[CANCER_GENES].sub(log_expr[REFERENCE].mean(axis=1), axis=0)
    return (x - x.min()) / (x.max() - x.min()) * 15


def score(s):
    """0–15 척도 발현 → 그룹 점수와 RS."""
    s = to_hgnc(s)
    out = pd.DataFrame({
        "her2_gs": (0.9 * s.GRB7 + 0.1 * s.ERBB2).clip(lower=8),
        "er_gs": (0.8 * s.ESR1 + 1.2 * s.PGR + s.BCL2 + s.SCUBE2) / 4,
        "proliferation_gs": s[PROLIFERATION].mean(axis=1).clip(lower=6.5),
        "invasion_gs": (s.CTSV + s.MMP11) / 2,
    })
    out["rs_unscaled"] = (0.47 * out.her2_gs - 0.34 * out.er_gs + 1.04 * out.proliferation_gs
                          + 0.10 * out.invasion_gs + 0.05 * s.CD68 - 0.08 * s.GSTM1 - 0.07 * s.BAG1)
    out["pseudo_rs"] = (20 * (out.rs_unscaled - 6.7)).clip(0, 100)
    return out


def assign_risk(rs):
    pct = rs.rank(pct=True, method="first")
    return pd.cut(pct, [0, LOW_PCT, 1 - HIGH_PCT, 1], labels=["low", "intermediate", "high"])


def build():
    m = pd.read_csv(MANIFEST).dropna(subset=["rna_file_id"])
    out = score(normalize(load_expression(m.patient_id)))
    out["rs_group"] = assign_risk(out.pseudo_rs)
    out["high_risk"] = out.rs_group == "high"
    out.index.name = "patient_id"
    out.round(4).to_csv(OUT)

    print("patients", len(out), "| pseudo_rs median", round(out.pseudo_rs.median(), 1))
    print(pd.crosstab(m.set_index("patient_id").group.reindex(out.index), out.rs_group, margins=True))
    print("saved", OUT)


if __name__ == "__main__":
    build()
