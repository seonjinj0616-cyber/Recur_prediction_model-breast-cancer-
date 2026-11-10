"""C2 (단계 5): P ∩ R 유전자로 원격 재발 점수를 METABRIC에서 학습하고 SCAN-B 검사 권고 코호트에서 검증한다.

1. 학습: METABRIC HR+/HER2− 재발 기록 완전 1,251명 (A4와 같은 환자). 유전자는 C1의 P ∩ R (c1_overlap.csv).
   발현은 코호트 안 유전자별 z-score (마이크로어레이 → RNA-seq 이동 대비). 사건 = DR, 10년에서 자름.
   유전자 수(243)가 사건 수(약 230)와 비슷해 L2(ridge) 규제 Cox. 규제 세기는 METABRIC 5-fold C-index로 고름.
   유전자를 고르지 않는 ridge를 쓰는 이유: 교집합 유전자는 모두 "H&E로 재현 가능한 재발 유전자"라 남기고,
   서로 강하게 얽힌 증식 유전자 묶음에서 대표 하나만 남기는 LASSO의 문제(A4)를 피함.
2. 고정: 가중치, 점수 0–100 변환(METABRIC 1–99 백분위), 위험 구간(METABRIC 점수 3등분), 10년 재발 확률(METABRIC 기준 생존).
3. 검증 (SCAN-B, 고정한 것을 바꾸지 않음): 10년 원격 재발 C-index와 부트스트랩 95% CI를 A3 · A5 기준선과 같은 환자로 비교,
   같은 개수의 무작위 유전자 세트 (N_RANDOM개, 같은 규제 세기), 위험 구간별 10년 KM 재발률,
   임상 저위험(NPI ≤ 3.4) 안에서 구간별 재발률, 예측 10년 확률 vs 관찰(KM) 5분위 보정.
실행: python -m module1_recurrence.models.scripts.c2_signature
결과: outputs/module1/c2_signature_weights.csv, c2_signature.json, data/{metabric,scanb}/c2_score.csv
"""
import json

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index
from sklearn.model_selection import StratifiedKFold

from module1_recurrence.data.scripts.expression import load, zscore
from module1_recurrence.data.scripts.external_fetch import DATA
from module1_recurrence.scripts.cohort import RECOMMENDED

OUT_DIR = DATA.parent / "outputs" / "module1"
TEN_YEARS = 3652.5
PENALIZERS = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0]
N_BOOT = 200
N_RANDOM = 100
NPI_LOW = 3.4  # Nottingham 예후 지수 "좋음" 구간 경계
GROUPS = ["low", "intermediate", "high"]


def survival_frame(cohort, expr):
    c = cohort.dropna(subset=["dr_event", "dr_days"]).set_index("patient_id")
    df = expr.reindex(c.index).fillna(0)  # z-score라 결측은 평균(0)
    df["T"] = c.dr_days.clip(upper=TEN_YEARS)
    df["E"] = ((c.dr_event == 1) & (c.dr_days <= TEN_YEARS)).astype(int)
    return df, c


def fit(df, genes, penalizer):
    return CoxPHFitter(penalizer=penalizer, l1_ratio=0.0).fit(df[[*genes, "T", "E"]], "T", "E")


def cindex(df, score):
    return concordance_index(df["T"], -np.asarray(score), df["E"])


def choose_penalizer(df, genes, seed=0):
    folds = list(StratifiedKFold(5, shuffle=True, random_state=seed).split(df, df["E"]))
    res = {}
    for pen in PENALIZERS:
        pred = pd.Series(index=df.index, dtype=float)
        for tr, te in folds:
            m = fit(df.iloc[tr], genes, pen)
            pred.iloc[te] = m.predict_log_partial_hazard(df.iloc[te][genes]).to_numpy()
        res[pen] = round(cindex(df, pred), 4)
    return max(res, key=res.get), res


def boot_ci(df, scores, rng):
    idx = [rng.integers(0, len(df), len(df)) for _ in range(N_BOOT)]
    out = {}
    for name, s in scores.items():
        s = np.asarray(s)
        boot = [concordance_index(df["T"].iloc[i], -s[i], df["E"].iloc[i]) for i in idx]
        out[name] = {"c_index": round(cindex(df, s), 3), "c_index_95ci": [round(float(x), 3) for x in np.percentile(boot, [2.5, 97.5])]}
    return out


def km_10y(df, mask):
    if mask.sum() == 0:
        return None
    k = KaplanMeierFitter().fit(df["T"][mask], df["E"][mask])
    return {"n": int(mask.sum()), "events": int(df["E"][mask].sum()),
            "recurrence_10y_%": round(100 * (1 - float(k.survival_function_at_times(TEN_YEARS - 1).iloc[0])), 1)}


def main():
    rng = np.random.default_rng(0)
    genes = pd.read_csv(OUT_DIR / "c1_overlap.csv").gene.tolist()

    mb_expr = zscore(load("metabric"))
    mb_cohort = pd.read_csv(DATA / "metabric" / "cohort.csv")
    mb_cohort = mb_cohort[mb_cohort.complete_rec]
    mb, mb_c = survival_frame(mb_cohort, mb_expr)
    genes = [g for g in genes if g in mb_expr.columns]

    pen, cv = choose_penalizer(mb, genes)
    model = fit(mb, genes, pen)
    raw_mb = model.predict_log_partial_hazard(mb[genes])
    lo, hi = np.percentile(raw_mb, [1, 99])
    to100 = lambda r: np.clip(100 * (np.asarray(r) - lo) / (hi - lo), 0, 100)  # noqa: E731
    cuts = np.percentile(raw_mb, [100 / 3, 200 / 3])
    groups = lambda r: pd.Series(np.digitize(np.asarray(r), cuts)).map(dict(enumerate(GROUPS))).to_numpy()  # noqa: E731
    prob10 = lambda d: 1 - model.predict_survival_function(d[genes], times=[TEN_YEARS - 1]).iloc[0].to_numpy()  # noqa: E731

    w = model.params_.rename("weight").to_frame().join(pd.read_csv(OUT_DIR / "c1_overlap.csv").set_index("gene")[["module", "direction", "r"]])
    w.sort_values("weight", key=abs, ascending=False).to_csv(OUT_DIR / "c2_signature_weights.csv")

    universe = pd.read_csv(OUT_DIR / "p_genes" / "p_genes_matched_dx1_tss.csv").gene
    universe = [g for g in universe if g in mb_expr.columns and g not in set(genes)]

    sb_expr = zscore(load("scanb"))
    sb_cohort = pd.read_csv(DATA / "scanb" / "cohort.csv")
    summary = {"genes": len(genes), "penalizer": pen, "metabric_cv_c_index": cv,
               "metabric_train": {"n": len(mb), "events_10y": int(mb["E"].sum())}}

    rand_sets = [list(rng.choice(universe, len(genes), replace=False)) for _ in range(N_RANDOM)]
    rand_models = [fit(mb, s, pen) for s in rand_sets]

    for name, expr, cohort in [("metabric", mb_expr, mb_cohort), ("scanb", sb_expr, sb_cohort)]:
        df, c = survival_frame(cohort, expr)
        raw = model.predict_log_partial_hazard(df[genes]).to_numpy()
        out = pd.DataFrame({"c2_raw": raw, "c2_score": to100(raw), "c2_group": groups(raw), "c2_prob10y": prob10(df)}, index=df.index)
        out.round(4).to_csv(DATA / name / "c2_score.csv")
        if name == "metabric":
            continue  # 학습 데이터: 점수만 저장 (교차검증 C-index는 위에서)

        base = pd.read_csv(DATA / name / "baselines.csv", index_col="patient_id")
        rec = c.group.isin(RECOMMENDED).to_numpy()
        res = {}
        for label, mask in [("recommended", rec), ("all HR+/HER2-", np.ones(len(df), bool))]:
            d = df[mask]
            s = out[mask].join(base)
            cols = ["c2_raw", "pseudo_rs", "proliferation", "sequoia", "npi", "ror"]
            ok = s[cols].notna().all(axis=1).to_numpy()
            d_ok, s_ok = d[ok], s[ok]
            r = {"n": int(ok.sum()), "events_10y": int(d_ok["E"].sum()), **boot_ci(d_ok, {k: s_ok[k] for k in cols}, rng)}
            rc = [cindex(d_ok, m.predict_log_partial_hazard(d_ok[st])) for m, st in zip(rand_models, rand_sets)]
            r["random_sets"] = {"mean": round(float(np.mean(rc)), 3), "95%": [round(float(x), 3) for x in np.percentile(rc, [2.5, 97.5])],
                                "c2_better_than_%": round(100 * float(np.mean(np.array(rc) < r["c2_raw"]["c_index"])), 1)}
            r["groups"] = {g: km_10y(d, (s.c2_group == g).to_numpy()) for g in GROUPS}
            npi_low = (s.npi <= NPI_LOW).to_numpy()
            r["npi_low_groups"] = {g: km_10y(d, npi_low & (s.c2_group == g).to_numpy()) for g in GROUPS}
            q = pd.qcut(s.c2_prob10y, 5, labels=False).to_numpy()
            r["calibration_quintiles"] = [{"predicted_%": round(100 * float(s.c2_prob10y[q == i].mean()), 1),
                                           "observed_%": km_10y(d, q == i)["recurrence_10y_%"]} for i in range(5)]
            res[label] = r
        summary[name] = res

    (OUT_DIR / "c2_signature.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
