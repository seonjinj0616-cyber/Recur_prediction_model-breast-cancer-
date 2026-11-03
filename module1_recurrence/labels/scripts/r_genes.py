"""A4 (단계 3): METABRIC에서 원격 재발과 연관된 유전자(R 목록)를 고른다.

환자: METABRIC HR+/HER2− 재발 기록 완전 1,251명. 사건 = DR, 시간 = TDR (전체 추적).
1) 유전자별 Cox: 유전자 z-score + 임상 보정(log1p 림프절 수, 종양 크기, 등급, 나이) → 표준편차당 위험비, BH FDR.
2) FDR < 0.05 유전자로 Cox elastic net(L1:L2 = 0.5:0.5, 임상 변수는 벌점 없음).
   LASSO만 쓰면 서로 연관된 증식 유전자 묶음에서 매번 하나만 골라 묶음 전체가 70%에 못 미친다. 벌점 크기는 5-fold CV C-index로 한 번 정하고 고정,
   부트스트랩 100회에서 선택된 비율 ≥ 0.7인 유전자를 R 목록으로 둔다.
결측: 유전자는 중앙값, 임상 변수도 중앙값으로 채운다.
실행: python -m module1_recurrence.labels.scripts.r_genes
결과: outputs/module1/r_gene_cox_all.csv (전체 유전자), outputs/module1/r_genes.csv (FDR 통과 유전자 + 선택 빈도)
"""
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from lifelines import CoxPHFitter
from scipy.stats import false_discovery_control
from sklearn.model_selection import KFold
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv

from module1_recurrence.data.scripts.expression import load, zscore
from module1_recurrence.data.scripts.external_fetch import DATA
from module1_recurrence.labels.scripts.pseudo_rs import CANCER_GENES, PROLIFERATION

OUT_DIR = DATA.parent / "outputs" / "module1"
L1_RATIO = 0.5
CLINICAL = ["log_nodes", "tumor_size_mm", "grade", "age"]
FDR = 0.05
N_BOOT, MIN_FREQ = 100, 0.7


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
    return out.sort_values("p")


def coxnet(X, y, n_clin, alphas=None):
    penalty = np.r_[np.zeros(n_clin), np.ones(X.shape[1] - n_clin)]
    return CoxnetSurvivalAnalysis(l1_ratio=L1_RATIO, penalty_factor=penalty, alphas=alphas,
                                  alpha_min_ratio=0.01, n_alphas=50, max_iter=100_000).fit(X, y)


def choose_alpha(X, y, n_clin):
    alphas = coxnet(X, y, n_clin).alphas_
    scores = np.zeros((5, len(alphas)))
    for k, (tr, te) in enumerate(KFold(5, shuffle=True, random_state=0).split(X)):
        m = coxnet(X[tr], y[tr], n_clin, alphas)
        scores[k] = [concordance_index_censored(y[te]["e"], y[te]["t"], m.predict(X[te], alpha=a))[0]
                     for a in alphas]
    best = int(scores.mean(axis=0).argmax())
    return alphas, best, scores.mean(axis=0)


def bootstrap_freq(X, y, n_clin, alphas, best, genes):
    rng = np.random.default_rng(0)
    path = alphas[:best + 1]  # 큰 벌점부터 따라 내려와야 수렴이 안정적이다

    def one(seed):
        i = np.random.default_rng(seed).integers(0, len(X), len(X))
        return coxnet(X[i], y[i], n_clin, path).coef_[n_clin:, -1] != 0

    picks = Parallel(n_jobs=-1)(delayed(one)(s) for s in rng.integers(0, 2**31, N_BOOT))
    return pd.Series(np.mean(picks, axis=0), index=genes)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    clin, expr, t, e = load_data()
    print(f"patients {len(t)}, distant recurrence {e.sum()}, genes {expr.shape[1]:,}")

    cox = screen(clin, expr, t, e)
    cox.to_csv(OUT_DIR / "r_gene_cox_all.csv")
    cand = cox[cox.fdr < FDR]
    print(f"유전자별 Cox FDR < {FDR}: {len(cand):,}개 (위험 {int((cand.hr_per_sd > 1).sum())} / 보호 {int((cand.hr_per_sd < 1).sum())})")

    X = pd.concat([clin, expr[cand.index]], axis=1).to_numpy()
    y = Surv.from_arrays(e.to_numpy().astype(bool), t.to_numpy(), "e", "t")
    alphas, best, cv = choose_alpha(X, y, len(CLINICAL))
    print(f"elastic net alpha {alphas[best]:.4g} (CV C-index {cv[best]:.3f})")

    cand = cand.assign(selection_freq=bootstrap_freq(X, y, len(CLINICAL), alphas, best, cand.index))
    cand["selected"] = cand.selection_freq >= MIN_FREQ
    cand["oncotype"] = cand.index.isin(CANCER_GENES)
    cand["oncotype_proliferation"] = cand.index.isin(PROLIFERATION)
    cand = cand.sort_values(["selection_freq", "p"], ascending=[False, True])
    cand.to_csv(OUT_DIR / "r_genes.csv")

    r = cand[cand.selected]
    print(f"R 목록 (선택 비율 ≥ {MIN_FREQ}): {len(r)}개, Oncotype 유전자 {int(r.oncotype.sum())}개 "
          f"(증식 {int(r.oncotype_proliferation.sum())}개)")
    print(r[["hr_per_sd", "fdr", "selection_freq"]].round(4).to_string())
    print(f"FDR 통과 Oncotype 유전자: {sorted(cand.index[cand.oncotype])}")
    print("saved", OUT_DIR / "r_genes.csv")


if __name__ == "__main__":
    main()
