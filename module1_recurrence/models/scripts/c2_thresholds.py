"""P 목록 r 기준별 재발 점수 비교: P(r ≥ 0.3 / 0.4 / 0.5) ∩ R로 C2를 각각 만들어 같은 방법으로 평가한다.

P 목록: 1,045명(HR+/HER2− + 다른 아형)으로 학습, HR+/HER2− 환자 안에서 잰 r (기관 단위 교차검증 예측) + SEQUOIA 기준 판정(전체 환자 기준).
기준선: 이전 513명 P 목록의 P ∩ R (243개, c1_overlap.csv).
각 묶음마다
1. METABRIC HR+/HER2− 1,251명 10년 원격 재발로 L2 Cox (c2_signature와 같은 방법, 규제 세기는 5-fold C-index로 선택) → 가중치 저장
2. SCAN-B 검사 권고 C-index (측정 RNA, 부트스트랩 95% CI)
3. TCGA HR+/HER2−: 유전자별 H&E 예측 발현(교차검증) → z-score → 같은 가중치 합 = H&E 점수.
   측정 RNA 점수와의 상관(재현도), 10년 원격 재발 C-index (H&E · 측정 RNA)
HistologyHSI-BC는 저장한 가중치로 서버 c2_hsi_scores를 돌려 따로 평가한다.
실행: python -m module1_recurrence.models.scripts.c2_thresholds
결과: outputs/module1/c2_thresholds/{name}_weights.csv, c2_thresholds.json
"""
import json

import numpy as np
import pandas as pd
from lifelines.utils import concordance_index

from module1_recurrence.data.scripts.expression import load, zscore
from module1_recurrence.data.scripts.external_fetch import DATA
from module1_recurrence.models.scripts.c2_signature import TEN_YEARS, choose_penalizer, fit, survival_frame
from module1_recurrence.scripts.cohort import RECOMMENDED

OUT = DATA.parent / "outputs" / "module1" / "c2_thresholds"
P_NEW = DATA.parent / "outputs" / "module1" / "p_genes_1045"
N_BOOT = 500


def hr_r(pred, expr, hr):
    y, p = expr.loc[hr, pred.columns], pred.loc[hr]
    yc, pc = y - y.mean(), p - p.mean()
    return (yc * pc).sum() / np.sqrt((yc ** 2).sum() * (pc ** 2).sum())


def ci(t, e, s, rng):
    t, e, s = map(np.asarray, (t, e, s))
    idx = [rng.integers(0, len(t), len(t)) for _ in range(N_BOOT)]
    boot = [concordance_index(t[i], -s[i], e[i]) for i in idx]
    return {"c_index": round(concordance_index(t, -s, e), 3), "95ci": [round(float(x), 3) for x in np.percentile(boot, [2.5, 97.5])]}


def main():
    rng = np.random.default_rng(0)
    OUT.mkdir(parents=True, exist_ok=True)
    new = pd.read_csv(P_NEW / "p_genes_matched_dx1_tss.csv").set_index("gene")
    pred = pd.read_parquet(P_NEW / "p_genes_matched_dx1_tss_pred.parquet")
    expr_all = pd.read_parquet(DATA / "tcga" / "all" / "expr.parquet")
    m_all = pd.read_csv(DATA / "tcga" / "all" / "manifest.csv").set_index("patient_id")
    hr = pred.index[pred.index.isin(m_all.index[m_all.cohort == "HR+/HER2-"])]
    r = hr_r(pred, expr_all, hr).reindex(new.index)
    R = set(pd.read_csv(OUT.parent / "r_genes.csv").gene)

    mb_expr, sb_expr = zscore(load("metabric")), zscore(load("scanb"))
    common = set(mb_expr.columns) & set(sb_expr.columns)
    sets = {"old_513_r0.3": pd.read_csv(OUT.parent / "c1_overlap.csv").gene.tolist()}
    for t in [0.3, 0.4, 0.5]:
        sets[f"hr1045_r{t}"] = sorted(set(new.index[new.well & (r >= t)]) & R)
    sets = {k: [g for g in v if g in common] for k, v in sets.items()}

    mbc = pd.read_csv(DATA / "metabric" / "cohort.csv")
    mb, _ = survival_frame(mbc[mbc.complete_rec], mb_expr)
    sb, sbc = survival_frame(pd.read_csv(DATA / "scanb" / "cohort.csv"), sb_expr)
    rec = sbc.group.isin(RECOMMENDED).to_numpy()

    # TCGA HR+/HER2−: 측정 RNA와 H&E 예측 발현 (둘 다 HR+ 환자 안 z-score), 원격 재발
    dist = pd.read_csv(DATA / "tcga" / "distant.csv", index_col=0).reindex(hr)
    tt = dist.dm_days.clip(upper=TEN_YEARS).to_numpy()
    te = ((dist.dm_event == 1) & (dist.dm_days <= TEN_YEARS)).astype(int).to_numpy()
    z = lambda d: ((d - d.mean()) / d.std()).replace([np.inf, -np.inf], np.nan).fillna(0)  # noqa: E731
    meas, he = z(expr_all.loc[hr]), z(pred.loc[hr])

    res = {"tcga_hr_patients": len(hr), "tcga_distant_10y": int(te.sum()), "scanb_recommended": int(rec.sum()),
           "scanb_events_10y": int(sb.E[rec].sum()), "sets": {}}
    for name, genes in sets.items():
        pen, cv = choose_penalizer(mb, genes)
        model = fit(mb, genes, pen)
        w = model.params_
        w.rename("weight").rename_axis("covariate").to_csv(OUT / f"{name}_weights.csv")
        sb_score = model.predict_log_partial_hazard(sb[genes]).to_numpy()[rec]
        g_he = [g for g in genes if g in he.columns]
        s_meas, s_he = meas[g_he] @ w[g_he], he[g_he] @ w[g_he]
        res["sets"][name] = {
            "genes": len(genes), "penalizer": pen, "metabric_cv_c_index": cv[pen],
            "scanb_recommended": ci(sb["T"][rec], sb.E[rec], sb_score, rng),
            "tcga_he_vs_rna_score_r": round(float(np.corrcoef(s_he, s_meas)[0, 1]), 3),
            "tcga_distant_he": ci(tt, te, s_he, rng), "tcga_distant_rna": ci(tt, te, s_meas, rng)}
        print(name, json.dumps(res["sets"][name], ensure_ascii=False), flush=True)
    (OUT / "c2_thresholds.json").write_text(json.dumps(res, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
