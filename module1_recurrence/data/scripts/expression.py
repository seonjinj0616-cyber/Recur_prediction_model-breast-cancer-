"""METABRIC · SCAN-B 발현 행렬을 코호트 환자만 남기고 환자 × 유전자 log2 값으로 정리한다.

METABRIC: Illumina HT-12 log2 강도 그대로, 같은 유전자 기호의 여러 프로브는 평균.
SCAN-B: 라이브러리 보정 FPKM → log2(FPKM + 0.1) (Staaf et al. 2022와 같은 변환), Ensembl ID → 유전자 기호.
세 데이터(METABRIC · SCAN-B · TCGA STAR-Counts)에 모두 있는 유전자만 남긴다.
z-score는 분석에 쓰는 환자 범위에 따라 달라지므로 저장하지 않고 zscore()로 그때 계산한다.

먼저 external_fetch --expression, external_cohort를 실행한다.
실행: python -m module1_recurrence.data.scripts.expression
결과: data/metabric/expr_hr_her2neg.parquet, data/scanb/expr_hr_her2neg.parquet (행: patient_id), data/expr_common_genes.txt
"""
import numpy as np
import pandas as pd

from module1_recurrence.data.scripts.external_fetch import DATA, METABRIC, SCANB

TCGA_RNA = DATA / "tcga" / "rnaseq"
OUT = {"metabric": DATA / "metabric" / "expr_hr_her2neg.parquet",
       "scanb": DATA / "scanb" / "expr_hr_her2neg.parquet"}
COMMON = DATA / "expr_common_genes.txt"


def tcga_genes():
    """TCGA STAR-Counts 파일의 유전자 이름 (GENCODE v36). 모든 파일이 같은 유전자 목록을 쓴다."""
    f = next(TCGA_RNA.glob("*.tsv"))
    d = pd.read_csv(f, sep="\t", comment="#", usecols=["gene_id", "gene_name"])
    return set(d[d.gene_id.str.startswith("ENSG")].gene_name)


def metabric_expr(cohort):
    m = pd.read_csv(METABRIC / "data_mrna_illumina_microarray.txt", sep="\t")
    m = m.drop(columns="Entrez_Gene_Id").groupby("Hugo_Symbol").mean()
    return m[cohort.expr_id].T.set_axis(cohort.patient_id, axis=0)


def scanb_expr(cohort):
    path = SCANB / "SCANB.9206.genematrix_noNeg.txt"
    gene_col = pd.read_csv(path, sep="\t", nrows=0).columns[0]  # 첫 열(Ensembl ID)은 이름이 비어 있다
    fpkm = pd.read_csv(path, sep="\t", index_col=0, usecols=[gene_col, *cohort.expr_id],
                       dtype={c: "float32" for c in cohort.expr_id})
    ann = pd.read_csv(SCANB / "Gene.ID.ann.txt", sep="\t", index_col="Gene.ID")["Gene.Name"]
    log = np.log2(fpkm[cohort.expr_id] + 0.1)
    log = log.groupby(ann.reindex(log.index).values).mean()
    return log.T.set_axis(cohort.patient_id, axis=0)


def zscore(expr):
    """유전자별로 주어진 환자 범위 안에서 평균 0, 표준편차 1."""
    return (expr - expr.mean()) / expr.std()


def load(name):
    """name: 'metabric' 또는 'scanb' → 환자 × 공통 유전자 log2 발현."""
    return pd.read_parquet(OUT[name])


def main():
    mb_cohort = pd.read_csv(DATA / "metabric" / "cohort.csv")
    sb_cohort = pd.read_csv(DATA / "scanb" / "cohort.csv")

    mb = metabric_expr(mb_cohort)
    sb = scanb_expr(sb_cohort)
    common = sorted(set(mb.columns) & set(sb.columns) & tcga_genes())
    print(f"genes: METABRIC {mb.shape[1]:,} | SCAN-B {sb.shape[1]:,} | 세 데이터 공통 {len(common):,}")

    for name, expr in [("metabric", mb), ("scanb", sb)]:
        expr = expr[common].astype("float32")
        expr.columns.name = None
        expr.to_parquet(OUT[name])
        print(f"{name}: {expr.shape[0]:,} patients x {expr.shape[1]:,} genes, NaN {int(expr.isna().sum().sum())} → {OUT[name]}")
    COMMON.write_text("\n".join(common), encoding="utf-8")


if __name__ == "__main__":
    main()
