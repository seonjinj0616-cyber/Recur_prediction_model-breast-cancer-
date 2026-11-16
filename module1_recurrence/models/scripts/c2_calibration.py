"""C2 재보정 (보강 3): METABRIC에서 만든 C2 점수의 10년 원격 재발 확률을 최근 코호트(SCAN-B)에 맞추고 위험 구간을 정한다.

문제: METABRIC(1980–2000년대 환자) 기준 생존으로 계산한 10년 확률이 SCAN-B 관찰보다 2.5–4배 높음 (C2 1차 결과).
방법: SCAN-B 검사 권고 코호트를 재발 여부 비율을 유지한 채 무작위로 반씩 나눔 (seed 0).
- 재보정 반: Cox(원격 재발 ~ C2 선형 점수) → 기울기와 기준 생존 S0(10년)을 다시 추정. 점수 계산식(유전자 가중치)은 그대로.
  10년 확률 = 1 − S0(10년) ^ exp(β × (점수 − 재보정 반 평균))
- 확인 반: 원래(METABRIC) 확률과 재보정 확률의 5분위별 예측 vs 관찰(KM), 전체 예측 평균 vs 관찰, 보정 기울기(Cox 계수),
  C-index(순위는 재보정으로 바뀌지 않음), 확률 기준 위험 구간(< 5%, 5–10%, ≥ 10%)별 관찰 10년 재발률.
위험 구간 기준은 탐색용 (Oncotype 저위험 · 고위험 구간의 10년 원격 재발이 대략 5% 미만 · 10% 이상인 것을 참고).
먼저 c2_signature를 실행해 data/scanb/c2_score.csv를 만든다.
실행: python -m module1_recurrence.models.scripts.c2_calibration
결과: outputs/module1/c2_calibration.json
"""
import json

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index
from sklearn.model_selection import train_test_split

from module1_recurrence.data.scripts.external_fetch import DATA
from module1_recurrence.scripts.cohort import RECOMMENDED

OUT = DATA.parent / "outputs" / "module1" / "c2_calibration.json"
TEN_YEARS = 3652.5
CUTS = [0.05, 0.10]
GROUPS = ["low (<5%)", "intermediate (5-10%)", "high (>=10%)"]


def km10(t, e):
    if len(t) == 0:
        return None
    k = KaplanMeierFitter().fit(t, e)
    return round(100 * (1 - float(k.survival_function_at_times(TEN_YEARS - 1).iloc[0])), 1)


def calib_table(d, prob):
    q = pd.qcut(prob, 5, labels=False)
    t, e = d["T"].to_numpy(), d.E.to_numpy()
    return [{"predicted_%": round(100 * float(prob[q == i].mean()), 1), "observed_%": km10(t[q == i], e[q == i])} for i in range(5)]


def main():
    c = pd.read_csv(DATA / "scanb" / "cohort.csv").set_index("patient_id")
    s = pd.read_csv(DATA / "scanb" / "c2_score.csv", index_col="patient_id")
    d = c[c.group.isin(RECOMMENDED)].dropna(subset=["dr_event", "dr_days"]).join(s, how="inner")
    d["T"] = d.dr_days.clip(upper=TEN_YEARS)
    d["E"] = ((d.dr_event == 1) & (d.dr_days <= TEN_YEARS)).astype(int)
    cal, val = train_test_split(d, test_size=0.5, stratify=d.E, random_state=0)

    mean = cal.c2_raw.mean()
    cox = CoxPHFitter().fit(pd.DataFrame({"x": cal.c2_raw - mean, "T": cal["T"], "E": cal.E}), "T", "E")
    beta = float(cox.params_["x"])
    s0 = float(cox.baseline_survival_.loc[:TEN_YEARS - 1].iloc[-1, 0])
    recal = lambda x: 1 - s0 ** np.exp(beta * (x - mean))  # noqa: E731

    p_orig, p_new = val.c2_prob10y.to_numpy(), recal(val.c2_raw.to_numpy())
    slope = CoxPHFitter().fit(pd.DataFrame({"lp": beta * (val.c2_raw - mean), "T": val["T"], "E": val.E}), "T", "E").params_["lp"]
    groups = pd.cut(p_new, [-np.inf, *CUTS, np.inf], labels=GROUPS, right=False)
    res = {"n_calibration": len(cal), "n_validation": len(val), "events_validation": int(val.E.sum()),
           "recalibration": {"beta_per_raw_unit": round(beta, 3), "S0_10y": round(s0, 4), "mean_raw": round(float(mean), 4)},
           "validation": {
               "observed_10y_%": km10(val["T"], val.E),
               "mean_predicted_original_%": round(100 * float(p_orig.mean()), 1),
               "mean_predicted_recalibrated_%": round(100 * float(p_new.mean()), 1),
               "calibration_slope": round(float(slope), 3),
               "c_index": round(concordance_index(val["T"], -val.c2_raw, val.E), 3),
               "quintiles_original": calib_table(val, p_orig),
               "quintiles_recalibrated": calib_table(val, p_new),
               "risk_groups": {g: {"n": int((groups == g).sum()), "events": int(val.E[groups == g].sum()),
                                   "observed_10y_%": km10(val["T"][groups == g], val.E[groups == g])} for g in GROUPS}}}
    OUT.write_text(json.dumps(res, indent=2, ensure_ascii=False))
    print(json.dumps(res, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
