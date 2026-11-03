"""A4 (단계 3): METABRIC에서 원격 재발과 연관된 유전자(R 목록)를 고른다.

환자: METABRIC HR+/HER2− 재발 기록 완전 1,251명. 사건 = DR, 시간 = TDR (전체 추적).
유전자별 Cox: 유전자 z-score + 임상 보정(log1p 림프절 수, 종양 크기, 등급, 나이) → 표준편차당 위험비, BH FDR.
R 목록 = 보정 FDR < 0.05인 유전자 전부. 단계 4에서 P 목록과 교집합을 내므로 연관된 유전자를 모두 남긴다.
(다변수 선택(LASSO · elastic net)은 서로 연관된 증식 유전자 묶음에서 대표만 남겨 이 목적에 맞지 않았다.)
비교용으로 임상 보정 없는 Cox 결과도 같이 저장한다.
결측: 유전자는 중앙값, 임상 변수도 중앙값으로 채운다.
실행: python -m module1_recurrence.labels.scripts.r_genes
결과: outputs/module1/r_gene_cox_all.csv (전체 유전자), outputs/module1/r_genes.csv (R 목록)
"""
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from lifelines import CoxPHFitter
from scipy.stats import false_discovery_control

from module1_recurrence.data.scripts.expression import load, zscore
from module1_recurrence.data.scripts.external_fetch import DATA
from module1_recurrence.labels.scripts.pseudo_rs import CANCER_GENES, PROLIFERATION

OUT_DIR = DATA.parent / "outputs" / "module1"
CLINICAL = ["log_nodes", "tumor_size_mm", "grade", "age"]
FDR = 0.05


def load_data():
    c = pd.read_csv(DATA / "metabric" / "cohort.csv")
    c = c[c.complete_rec].set_index("patient_id")
    clin = c.assign(log_nodes=np.log1p(c.nodes))[CLINICAL]
    clin = zscore(clin.fillna(clin.median()))
    expr = load("metabric").loc[c.index]
    expr = zscore(expr.fillna(expr.median()))
    return clin, expr, c.dr_days, c.dr_event.astype(int)


def gene_cox(gene, x, clin, t, e):
    d = clin.assign(gene=x, t=t, e=e)
    s = CoxPHFitter().fit(d, "t", "e").summary.loc["gene"]
    return gene, s["exp(coef)"], s["exp(coef) lower 95%"], s["exp(coef) upper 95%"], s.p


def screen(clin, expr, t, e):
    rows = Parallel(n_jobs=-1, batch_size=64)(delayed(gene_cox)(g, expr[g], clin, t, e) for g in expr.columns)
    out = pd.DataFrame(rows, columns=["gene", "hr_per_sd", "hr_lower", "hr_upper", "p"]).set_index("gene")
    out["fdr"] = false_discovery_control(out.p)
    return out


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    clin, expr, t, e = load_data()
    print(f"patients {len(t)}, distant recurrence {e.sum()}, genes {expr.shape[1]:,}")

    cox = screen(clin, expr, t, e)
    unadj = screen(clin[[]], expr, t, e)[["hr_per_sd", "fdr"]].add_suffix("_unadj")
    cox = cox.join(unadj).sort_values("p")
    cox["oncotype"] = cox.index.isin(CANCER_GENES)
    cox["oncotype_proliferation"] = cox.index.isin(PROLIFERATION)
    cox.to_csv(OUT_DIR / "r_gene_cox_all.csv")

    r = cox[cox.fdr < FDR]
    r.to_csv(OUT_DIR / "r_genes.csv")
    lost = int(((cox.fdr_unadj < FDR) & (cox.fdr >= FDR)).sum())
    print(f"보정 없음 FDR < {FDR}: {int((cox.fdr_unadj < FDR).sum()):,}개 → 임상 보정 후 {len(r)}개 "
          f"(보정으로 빠짐 {lost:,}개)")
    print(f"R 목록 {len(r)}개: 위험 {int((r.hr_per_sd > 1).sum())} / 보호 {int((r.hr_per_sd < 1).sum())}, "
          f"Oncotype 유전자 {sorted(r.index[r.oncotype])}")
    print("saved", OUT_DIR / "r_genes.csv")


if __name__ == "__main__":
    main()
