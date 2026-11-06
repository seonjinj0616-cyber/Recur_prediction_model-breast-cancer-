"""TCGA 코호트(manifest 515명)의 RNA-seq을 환자 × 공통 유전자 log2(TPM + 1) 행렬로 정리한다 (B4 정답).

STAR-Counts의 tpm_unstranded를 쓰고, 같은 유전자 기호가 여러 Ensembl ID에 걸리면 더한다 (pseudo_rs와 같음).
유전자는 METABRIC · SCAN-B · TCGA 공통 17,812개 (data/expr_common_genes.txt, expression.py가 만든 것).
먼저 rnaseq_download, expression을 실행한다.
실행: python -m module1_recurrence.data.scripts.tcga_expression
결과: data/tcga/expr.parquet (행: patient_id)
"""
import numpy as np
import pandas as pd

from module1_recurrence.data.scripts.expression import COMMON, TCGA_RNA

OUT = TCGA_RNA.parent / "expr.parquet"


def main():
    genes = COMMON.read_text().split()
    m = pd.read_csv(TCGA_RNA.parent / "manifest.csv")
    cols = {}
    for pid in m.patient_id:
        d = pd.read_csv(TCGA_RNA / f"{pid}.tsv", sep="\t", comment="#", usecols=["gene_name", "tpm_unstranded"])
        cols[pid] = d.groupby("gene_name").tpm_unstranded.sum().reindex(genes)
    expr = np.log2(pd.DataFrame(cols).T.astype("float32") + 1)
    expr.index.name = "patient_id"
    missing = expr.columns[expr.isna().all()].tolist()
    expr.to_parquet(OUT)
    print(expr.shape, "→", OUT, "| 모든 환자에서 값 없는 유전자:", missing[:10], len(missing))


if __name__ == "__main__":
    main()
