"""C1 (단계 4): P 목록(H&E로 잘 예측) ∩ R 목록(METABRIC 원격 재발 연관)과 그 구성.

1. P = B4 결과에서 SEQUOIA 기준 "잘 예측됨" + r ≥ 0.3. 공통 유전자(B4 대상)만 놓고 R과 교집합.
   우연 기대값과 초기하 검정 (모집단 = B4 대상 유전자 중 R 판정이 가능한 것).
2. 유전자 묶음: TCGA 발현에서 표지 유전자 묶음 점수(묶음 안 z-score 평균)와의 Pearson r이 가장 큰 묶음,
   단 |r| ≥ 0.5일 때만 배정하고 나머지는 "기타". 표지 유전자 자신은 그 묶음으로 둔다.
   묶음: 증식 · 면역 · 기질 · 에스트로겐 반응 (표지는 아래 MODULES, 널리 쓰이는 대표 유전자 5개씩).
3. R의 방향(위험비 > 1: 재발 위험 ↑ / < 1: ↓)과 묶음별 개수를 P ∩ R · R 전체 · P 전체 · 전체 유전자에서 비교.
실행: python -m module1_recurrence.models.scripts.c1_overlap [--p outputs/module1/p_genes/p_genes_matched_dx1_tss.csv]
결과: outputs/module1/c1_overlap.csv (교집합 유전자), c1_overlap.json (요약)
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from module1_recurrence.data.scripts.external_fetch import DATA

OUT_DIR = DATA.parent / "outputs" / "module1"
MODULES = {
    "proliferation": ["MKI67", "TOP2A", "BUB1", "CCNB1", "AURKA"],
    "immune": ["PTPRC", "CD3E", "CD2", "CD52", "LCP2"],
    "stroma": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"],
    "estrogen": ["ESR1", "PGR", "GATA3", "FOXA1", "XBP1"],
}
MIN_R = 0.3
MIN_MODULE_R = 0.5


def module_of(expr):
    """유전자 → 가장 가까운 표지 묶음 (|r| ≥ MIN_MODULE_R), 아니면 'other'."""
    z = (expr - expr.mean()) / expr.std()
    scores = pd.DataFrame({k: z[g].mean(1) for k, g in MODULES.items()})
    zc = z - z.mean()
    sc = (scores - scores.mean()) / scores.std(ddof=0)
    r = (zc.T @ sc) / (len(z) * zc.std(ddof=0).to_numpy()[:, None])
    best = r.abs().idxmax(axis=1).where(r.abs().max(axis=1) >= MIN_MODULE_R, "other")
    for k, genes in MODULES.items():
        best[genes] = k
    return best, r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=Path, default=OUT_DIR / "p_genes" / "p_genes_matched_dx1_tss.csv", help="B4 결과 csv")
    ap.add_argument("--r", type=Path, default=OUT_DIR / "r_genes.csv")
    ap.add_argument("--r_all", type=Path, default=OUT_DIR / "r_gene_cox_all.csv", help="R 판정이 가능한 전체 유전자")
    ap.add_argument("--expr", type=Path, default=DATA / "tcga" / "expr.parquet")
    args = ap.parse_args()

    p = pd.read_csv(args.p).set_index("gene")
    r = pd.read_csv(args.r).set_index("gene")
    universe = p.index.intersection(pd.read_csv(args.r_all).gene)
    p = p.loc[universe]
    P = set(p.index[p.well & (p.r >= MIN_R)])
    R = set(r.index) & set(universe)
    both = sorted(P & R)
    expected = len(P) * len(R) / len(universe)
    pval = hypergeom.sf(len(both) - 1, len(universe), len(R), len(P))

    expr = pd.read_parquet(args.expr)[universe]
    module, _ = module_of(expr)
    sets = {"P∩R": both, "R": sorted(R), "P": sorted(P), "all": list(universe)}
    comp = {k: module[g].value_counts().reindex([*MODULES, "other"], fill_value=0).to_dict() for k, g in sets.items()}

    out = r.loc[both, ["hr_per_sd", "fdr"]].join(p.loc[both, ["rank", "r"]]).assign(module=module[both])
    out["direction"] = np.where(out.hr_per_sd > 1, "risk", "protective")
    out = out.sort_values("r", ascending=False)
    out.to_csv(OUT_DIR / "c1_overlap.csv")

    summary = {"p_file": args.p.name, "universe": len(universe), "P": len(P), "R": len(R), "P∩R": len(both),
               "expected": round(expected, 1), "fold": round(len(both) / expected, 2), "hypergeom_p": float(f"{pval:.2e}"),
               "R_in_P_frac": round(len(both) / len(R), 3), "direction_P∩R": out.direction.value_counts().to_dict(),
               "direction_R": pd.Series(np.where(r.loc[sorted(R)].hr_per_sd > 1, "risk", "protective")).value_counts().to_dict(),
               "modules": comp}
    (OUT_DIR / "c1_overlap.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(out.groupby(["module", "direction"]).size().unstack(fill_value=0))
    print(out.head(15).round(3).to_string())


if __name__ == "__main__":
    main()
