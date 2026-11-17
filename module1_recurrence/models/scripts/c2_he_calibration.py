"""H&E 점수 전용 재보정: H&E로 예측한 C2 점수를 10년 원격 재발 확률 · 위험 구간으로 바꾸는 방법을 비교한다.

배경: c2_calibration은 측정 RNA 점수 → 확률 환산(SCAN-B 재보정 반의 Cox 기울기 β · 기준 생존 S0)을 만들었다.
H&E 점수는 Ridge 예측이라 RNA 점수보다 퍼짐이 작아 같은 환산을 쓰면 확률이 가운데로 몰린다.
환산은 단조 변환이라 순위(C-index)는 어느 방법이든 같고, 바뀌는 것은 확률값 · 위험 구간 배정이다.
방법
- naive: RNA 환산을 H&E 점수에 그대로 적용
- stretch: H&E 점수를 RNA 점수와 같은 평균 · 표준편차로 늘린 뒤 적용 (TCGA에서 맞춤)
- conditional: TCGA에서 RNA 점수 | H&E 점수 ~ N(a + b·H&E, σ²)을 추정하고, 10년 확률 = 1 − E[S0^exp(β(RNA − m))]
  (가우스-에르미트 적분). H&E로 RNA 점수를 정확히 알 수 없다는 불확실성을 확률에 반영하는 방법.
평가
1. 시뮬레이션 (SCAN-B 확인 반, 원격 재발 있음): TCGA에서 잰 H&E | RNA 관계(H&E = c + d·RNA + 잡음)로 각 환자의
   "H&E 같은 점수"를 만들고(10회 반복 평균), 방법별 예측 평균 · 5분위 예측 vs 관찰 · 위험 구간별 관찰 재발률 · C-index.
2. TCGA HR+/HER2− 실제: 방법별 H&E 위험 구간이 측정 RNA 위험 구간과 일치하는 비율, 구간별 환자 비율.
3. HistologyHSI-BC: 방법별 위험 구간 분포 (환자-대조 설계라 절대 확률은 평가하지 않음).
먼저 c2_signature, c2_calibration, c2_hsi_scores · c2_hsi_eval(data/hsi/c2_validation_scores.csv)을 실행한다.
실행: python -m module1_recurrence.models.scripts.c2_he_calibration
결과: outputs/module1/c2_he_calibration.json
"""
import json

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index
from sklearn.model_selection import train_test_split

from module1_recurrence.data.scripts.external_fetch import DATA
from module1_recurrence.scripts.cohort import RECOMMENDED

OUT = DATA.parent / "outputs" / "module1" / "c2_he_calibration.json"
TEN_YEARS = 3652.5
CUTS = [0.05, 0.10]
GROUPS = ["low (<5%)", "intermediate (5-10%)", "high (>=10%)"]
N_SIM = 10
GH_X, GH_W = np.polynomial.hermite_e.hermegauss(40)  # 표준정규 적분 노드


def km10(t, e):
    if len(t) == 0 or e.sum() == 0 and len(t) < 5:
        return None
    k = KaplanMeierFitter().fit(t, e)
    return round(100 * (1 - float(k.survival_function_at_times(TEN_YEARS - 1).iloc[0])), 1)


def group(p):
    return pd.cut(p, [-np.inf, *CUTS, np.inf], labels=GROUPS, right=False)


def scanb_recalibration():
    """c2_calibration과 같은 분할(seed 0)로 재보정 반에서 β · S0 · 평균, 확인 반 반환."""
    c = pd.read_csv(DATA / "scanb" / "cohort.csv").set_index("patient_id")
    s = pd.read_csv(DATA / "scanb" / "c2_score.csv", index_col="patient_id")
    d = c[c.group.isin(RECOMMENDED)].dropna(subset=["dr_event", "dr_days"]).join(s, how="inner")
    d["T"] = d.dr_days.clip(upper=TEN_YEARS)
    d["E"] = ((d.dr_event == 1) & (d.dr_days <= TEN_YEARS)).astype(int)
    cal, val = train_test_split(d, test_size=0.5, stratify=d.E, random_state=0)
    mean = cal.c2_raw.mean()
    cox = CoxPHFitter().fit(pd.DataFrame({"x": cal.c2_raw - mean, "T": cal["T"], "E": cal.E}), "T", "E")
    return float(cox.params_["x"]), float(cox.baseline_survival_.loc[:TEN_YEARS - 1].iloc[-1, 0]), float(mean), val


def main():
    rng = np.random.default_rng(0)
    beta, s0, mean, val = scanb_recalibration()
    rna_prob = lambda x: 1 - s0 ** np.exp(beta * (np.asarray(x) - mean))  # noqa: E731

    s = pd.read_csv(DATA / "hsi" / "c2_validation_scores.csv", index_col="id", dtype={"id": str})
    tc = s[s.cohort == "TCGA HR+/HER2- (CV)"]
    he, rna = tc.score.to_numpy(), tc.measured.to_numpy()
    # TCGA는 측정 RNA 점수를 TCGA 안 z-score로 만들었으므로 SCAN-B 점수와 같은 척도(Σ 가중치 × z)로 봄
    b, a = np.polyfit(he, rna, 1)
    sigma = float(np.std(rna - (a + b * he), ddof=2))
    d_, c_ = np.polyfit(rna, he, 1)
    s_noise = float(np.std(he - (c_ + d_ * rna), ddof=2))
    mu_he, sd_he, mu_rna, sd_rna = he.mean(), he.std(), rna.mean(), rna.std()

    methods = {
        "naive": lambda h: rna_prob(h),
        "stretch": lambda h: rna_prob(mu_rna + (np.asarray(h) - mu_he) * sd_rna / sd_he),
        "conditional": lambda h: 1 - (GH_W[None, :] / GH_W.sum() * s0 ** np.exp(
            beta * ((a + b * np.asarray(h))[:, None] + sigma * GH_X[None, :] - mean))).sum(1),
    }
    res = {"scanb_recalibration": {"beta": round(beta, 3), "S0_10y": round(s0, 4)},
           "tcga_relation": {"rna_given_he": {"a": round(a, 3), "b": round(b, 3), "sigma": round(sigma, 3)},
                             "he_given_rna": {"c": round(c_, 3), "d": round(d_, 3), "noise_sd": round(s_noise, 3)},
                             "r": round(float(np.corrcoef(he, rna)[0, 1]), 3), "sd_he": round(float(sd_he), 3), "sd_rna": round(float(sd_rna), 3)}}

    # 1. 시뮬레이션: SCAN-B 확인 반에 H&E 같은 점수
    t, e = val["T"].to_numpy(), val.E.to_numpy()
    sim = {k: {"mean_pred": [], "quint": [], "groups": [], "c": []} for k in methods}
    for _ in range(N_SIM):
        h = c_ + d_ * val.c2_raw.to_numpy() + rng.normal(0, s_noise, len(val))
        q = pd.qcut(h, 5, labels=False)
        for k, f in methods.items():
            p = f(h)
            g = group(p)
            sim[k]["mean_pred"].append(100 * p.mean())
            sim[k]["quint"].append([(100 * p[q == i].mean(), km10(t[q == i], e[q == i])) for i in range(5)])
            sim[k]["groups"].append([((g == G).mean() * 100, km10(t[g == G], e[g == G])) for G in GROUPS])
            sim[k]["c"].append(concordance_index(t, -h, e))
    res["simulation_scanb_validation_half"] = {"observed_10y_%": km10(t, e), "rna_score_reference": {
        "mean_pred_%": round(100 * float(rna_prob(val.c2_raw).mean()), 1),
        "c_index": round(concordance_index(t, -val.c2_raw, e), 3)}}
    for k, v in sim.items():
        avg = lambda xs: None if any(x is None for x in xs) else round(float(np.mean(xs)), 1)  # noqa: E731
        res["simulation_scanb_validation_half"][k] = {
            "mean_pred_%": round(float(np.mean(v["mean_pred"])), 1), "c_index": round(float(np.mean(v["c"])), 3),
            "quintiles_pred_vs_obs_%": [[avg([r[i][0] for r in v["quint"]]), avg([r[i][1] for r in v["quint"]])] for i in range(5)],
            "groups_share_and_observed_%": {G: [avg([r[j][0] for r in v["groups"]]), avg([r[j][1] for r in v["groups"]])]
                                            for j, G in enumerate(GROUPS)}}

    # 2. TCGA 실제: H&E 구간 vs 측정 RNA 구간
    g_rna = group(rna_prob(rna))
    res["tcga_group_agreement"] = {"rna_groups_%": {G: round(100 * float((g_rna == G).mean()), 1) for G in GROUPS}}
    for k, f in methods.items():
        g_he = group(f(he))
        res["tcga_group_agreement"][k] = {"agreement_%": round(100 * float((g_he == g_rna).mean()), 1),
                                          "groups_%": {G: round(100 * float((g_he == G).mean()), 1) for G in GROUPS}}

    # 3. HistologyHSI-BC 구간 분포
    clin = pd.read_excel(DATA / "hsi" / "clinical.xlsx", sheet_name="Clinical Data")
    clin.index = clin["Case ID"].astype(str)
    h = s[s.cohort == "HistologyHSI-BC"].drop(columns=["time", "event"], errors="ignore").join(clin, how="inner")
    h = h[((h.ER == 1) | (h.PR == 1)) & (h.HER2 == 0)]
    res["hsi_hr_groups"] = {}
    for k, f in methods.items():
        g = group(f(h.score.to_numpy()))
        res["hsi_hr_groups"][k] = {G: {"n": int((g == G).sum()), "relapse": int(h.Relapse.to_numpy()[g == G].sum())} for G in GROUPS}
    OUT.write_text(json.dumps(res, indent=2, ensure_ascii=False))
    print(json.dumps(res, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
