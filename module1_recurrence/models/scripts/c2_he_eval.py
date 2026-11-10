"""C2 H&E 평가: 유전자 묶음별 재발 점수를 METABRIC에서 학습하고, TCGA에서 H&E로 예측한 발현으로 계산해 PFI와 비교한다.

묶음 (모두 B4 대상 12,758개 안): P ∩ R (C2) / R / P (H&E로 잘 예측, r ≥ 0.3) / 전체 12,758개.
학습: METABRIC HR+/HER2− 재발 기록 완전 1,251명, 10년 원격 재발, 거의 ridge인 Coxnet (l1_ratio 0.01),
      규제 세기는 METABRIC 5-fold C-index로 선택. 묶음만 바꾸고 방법은 같음.
평가:
- SCAN-B 검사 권고 (측정 RNA): 묶음 자체의 예후력
- TCGA 맞춤 방식 482명: ① 측정 RNA (log2 TPM) ② H&E 예측 발현 (p_genes --save_pred, 기관 단위 교차검증 예측) 로 점수 계산
  → 10년 PFI C-index, 부트스트랩 95% CI, 같은 표본으로 묶음 간 차이.
  발현은 각 데이터 안에서 유전자별 z-score (H&E 예측값은 분산이 줄어 있어 척도를 맞춤).
  TCGA PFI는 국소 · 원격 재발 · 진행을 합친 것이라 METABRIC 원격 재발과 정의가 다르고, 사건이 적어 CI가 넓음.
실행: python -m module1_recurrence.models.scripts.c2_he_eval
결과: outputs/module1/c2_he_eval.json
"""
import json

import numpy as np
import pandas as pd
from lifelines.utils import concordance_index
from sklearn.model_selection import StratifiedKFold
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.util import Surv

from module1_recurrence.data.scripts.expression import load, zscore
from module1_recurrence.data.scripts.external_fetch import DATA
from module1_recurrence.scripts.cohort import RECOMMENDED

OUT_DIR = DATA.parent / "outputs" / "module1"
P_DIR = OUT_DIR / "p_genes"
TEN_YEARS = 3652.5
L1 = 0.01
N_BOOT = 500
MIN_R = 0.3


def z(expr):
    return zscore(expr).replace([np.inf, -np.inf], np.nan).fillna(0)


def frame(cohort, expr, time, event):
    c = cohort.dropna(subset=[time, event]).set_index("patient_id")
    c = c[c.index.isin(expr.index)]
    t = c[time].clip(upper=TEN_YEARS).to_numpy()
    e = ((c[event] == 1) & (c[time] <= TEN_YEARS)).to_numpy()
    return expr.loc[c.index], t, e, c


def cindex(t, e, s):
    return concordance_index(t, -np.asarray(s), e)


def train(x, t, e, seed=0):
    X, y = x.to_numpy(), Surv.from_arrays(e, t)
    path = CoxnetSurvivalAnalysis(l1_ratio=L1, alpha_min_ratio=0.01, n_alphas=30, max_iter=200000).fit(X, y).alphas_
    cv = np.zeros(len(path))
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=seed).split(X, e):
        m = CoxnetSurvivalAnalysis(l1_ratio=L1, alphas=path, max_iter=200000).fit(X[tr], y[tr])
        cv += [cindex(t[te], e[te], m.predict(X[te], alpha=a)) for a in path]
    best = path[int(np.argmax(cv))]
    model = CoxnetSurvivalAnalysis(l1_ratio=L1, alphas=[best], max_iter=200000).fit(X, y)
    return model, round(float(cv.max() / 5), 3)


def evaluate(t, e, scores, rng):
    """같은 부트스트랩 표본으로 묶음마다 C-index · CI와 P ∩ R 대비 차이."""
    idx = [rng.integers(0, len(t), len(t)) for _ in range(N_BOOT)]
    boot = {k: np.array([cindex(t[i], e[i], s[i]) for i in idx]) for k, s in scores.items()}
    ref = "P∩R"
    out = {}
    for k, s in scores.items():
        d = boot[k] - boot[ref]
        out[k] = {"c_index": round(cindex(t, e, s), 3), "95ci": [round(float(x), 3) for x in np.percentile(boot[k], [2.5, 97.5])],
                  "minus_P∩R_95ci": [round(float(x), 3) for x in np.percentile(d, [2.5, 97.5])]}
    return out


def main():
    rng = np.random.default_rng(0)
    p = pd.read_csv(P_DIR / "p_genes_matched_dx1_tss.csv")
    pred = pd.read_parquet(P_DIR / "p_genes_matched_dx1_tss_pred.parquet")
    mb_raw, sb_raw = load("metabric"), load("scanb")
    universe = [g for g in p.gene if g in mb_raw.columns and g in sb_raw.columns and g in pred.columns]
    P = set(p.gene[p.well & (p.r >= MIN_R)])
    sets = {"P∩R": pd.read_csv(OUT_DIR / "c1_overlap.csv").gene, "R": pd.read_csv(OUT_DIR / "r_genes.csv").gene,
            "P": sorted(P), "all": universe}
    sets = {k: [g for g in v if g in set(universe)] for k, v in sets.items()}

    mbc = pd.read_csv(DATA / "metabric" / "cohort.csv")
    mbx, mbt, mbe, _ = frame(mbc[mbc.complete_rec], z(mb_raw[universe]), "dr_days", "dr_event")
    sbx, sbt, sbe, sbc = frame(pd.read_csv(DATA / "scanb" / "cohort.csv"), z(sb_raw[universe]), "dr_days", "dr_event")
    rec = sbc.group.isin(RECOMMENDED).to_numpy()

    tc = pd.read_csv(DATA / "tcga" / "manifest.csv")
    measured = pd.read_parquet(DATA / "tcga" / "expr.parquet")
    pids = pred.index.intersection(measured.index)
    tmx, tt, te, tcc = frame(tc, z(measured.loc[pids, universe]), "PFI_time", "PFI")
    thx = z(pred.loc[tmx.index, universe])
    trec = tcc.recommended.astype(bool).to_numpy() if "recommended" in tcc else tcc.group.isin(RECOMMENDED).to_numpy()

    summary = {"sets": {k: len(v) for k, v in sets.items()}, "metabric_cv": {},
               "scanb_recommended": {"n": int(rec.sum()), "events_10y": int(sbe[rec].sum())},
               "tcga": {"n": len(tt), "events_10y": int(te.sum()), "recommended_n": int(trec.sum()),
                        "recommended_events_10y": int(te[trec].sum())}}
    sb_s, tm_s, th_s = {}, {}, {}
    for k, genes in sets.items():
        model, cv = train(mbx[genes], mbt, mbe)
        summary["metabric_cv"][k] = cv
        sb_s[k] = model.predict(sbx[genes].to_numpy())[rec]
        tm_s[k] = model.predict(tmx[genes].to_numpy())
        th_s[k] = model.predict(thx[genes].to_numpy())
        print(k, len(genes), "genes, METABRIC cv", cv, flush=True)

    summary["scanb_recommended"].update(evaluate(sbt[rec], sbe[rec], sb_s, rng))
    summary["tcga_measured_rna"] = evaluate(tt, te, tm_s, rng)
    summary["tcga_he_predicted"] = evaluate(tt, te, th_s, rng)
    summary["tcga_he_predicted_recommended"] = evaluate(tt[trec], te[trec], {k: v[trec] for k, v in th_s.items()}, rng)
    summary["tcga_measured_vs_he_score_r"] = {k: round(float(np.corrcoef(tm_s[k], th_s[k])[0, 1]), 3) for k in sets}
    (OUT_DIR / "c2_he_eval.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
