"""C2 + 임상 결합: 재발 점수에 임상 정보를 더하면 원격 재발 예측이 얼마나 좋아지는지 본다.

결합 가중치는 METABRIC HR+/HER2− 1,251명(10년 원격 재발)에서 Cox로 학습해 고정하고, 평가 데이터에서는 다시 맞추지 않는다.
점수는 데이터마다 표준편차 1로 맞춘다 (METABRIC: 측정 RNA C2, 평가: RNA 또는 H&E 점수).
- NPI 결합 (SCAN-B, HistologyHSI-BC): Cox(점수 + NPI). NPI = 0.2 × 크기(cm) + 등급 + 림프절 단계.
- 임상 결합 (TCGA, 등급이 없어 NPI 대신): Cox(점수 + T 병기 + log(1 + 양성 림프절) + 나이) — c2_clinical과 같은 변수.
평가: 임상만 · 점수만 · 점수 + 임상의 10년 원격 재발 C-index와, 같은 부트스트랩 표본에서 (점수 + 임상) − (임상만) 차이의 95% CI.
- SCAN-B 검사 권고 (측정 RNA C2)
- HistologyHSI-BC HR+/HER2− (H&E 점수, 주 평가) · 47명 전체(보조)
- TCGA HR+/HER2− (H&E 교차검증 점수)
먼저 c2_signature(c2_score.csv), c2_hsi_scores · c2_hsi_eval(data/hsi/c2_validation_scores.csv)을 실행한다.
실행: python -m module1_recurrence.models.scripts.c2_combined
결과: outputs/module1/c2_combined.json
"""
import json

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

from module1_recurrence.data.scripts.external_fetch import DATA
from module1_recurrence.models.scripts.c2_clinical import CLIN, std, t_from_size, tcga_clinical
from module1_recurrence.scripts.cohort import RECOMMENDED

OUT = DATA.parent / "outputs" / "module1" / "c2_combined.json"
TEN_YEARS = 3652.5
N_BOOT = 1000


def metabric_models():
    mb = pd.read_csv(DATA / "metabric" / "cohort.csv")
    mb = mb[mb.complete_rec].set_index("patient_id").join(pd.read_csv(DATA / "metabric" / "c2_score.csv", index_col="patient_id").c2_raw)
    d = pd.DataFrame({"score": std(mb.c2_raw), "npi": mb.npi, "t_stage": t_from_size(mb.tumor_size_mm), "log_nodes": np.log1p(mb.nodes),
                      "age": mb.age, "T": mb.dr_days.clip(upper=TEN_YEARS),
                      "E": ((mb.dr_event == 1) & (mb.dr_days <= TEN_YEARS)).astype(int)})
    d = d.fillna(d.median())
    fit = lambda cols: CoxPHFitter().fit(d[[*cols, "T", "E"]], "T", "E")  # noqa: E731
    return {"npi": fit(["npi"]), "score+npi": fit(["score", "npi"]), "clinical": fit(CLIN), "score+clinical": fit(["score", *CLIN])}


def evaluate(t, e, scores, base, rng):
    t, e = np.asarray(t), np.asarray(e)
    idx = [rng.integers(0, len(t), len(t)) for _ in range(N_BOOT)]
    idx = [i for i in idx if e[i].sum() > 0]
    out = {"n": int(len(t)), "events": int(e.sum())}
    boot = {}
    for k, s in scores.items():
        s = np.asarray(s)
        boot[k] = np.array([concordance_index(t[i], -s[i], e[i]) for i in idx])
        out[k] = {"c_index": round(concordance_index(t, -s, e), 3), "95ci": [round(float(x), 3) for x in np.percentile(boot[k], [2.5, 97.5])]}
    for k in scores:
        if "+" in k:
            d = boot[k] - boot[base]
            out[f"{k} minus {base}"] = {"mean": round(float(d.mean()), 3), "95ci": [round(float(x), 3) for x in np.percentile(d, [2.5, 97.5])]}
    return out


def main():
    rng = np.random.default_rng(0)
    mods = metabric_models()
    res = {"metabric_hazard_ratios": {k: m.hazard_ratios_.round(3).to_dict() for k, m in mods.items()}}

    # SCAN-B 검사 권고 (측정 RNA 점수)
    c = pd.read_csv(DATA / "scanb" / "cohort.csv").set_index("patient_id")
    d = c[c.group.isin(RECOMMENDED)].dropna(subset=["dr_event", "dr_days"]).join(
        pd.read_csv(DATA / "scanb" / "c2_score.csv", index_col="patient_id").c2_raw).join(
        pd.read_csv(DATA / "scanb" / "baselines.csv", index_col="patient_id").npi).dropna(subset=["c2_raw", "npi"])
    x = pd.DataFrame({"score": std(d.c2_raw), "npi": d.npi})
    res["SCAN-B recommended (RNA)"] = evaluate(
        d.dr_days.clip(upper=TEN_YEARS), ((d.dr_event == 1) & (d.dr_days <= TEN_YEARS)).astype(int),
        {"npi": d.npi, "score": d.c2_raw, "score+npi": mods["score+npi"].predict_log_partial_hazard(x)}, "npi", rng)

    # HistologyHSI-BC (H&E 점수)
    s = pd.read_csv(DATA / "hsi" / "c2_validation_scores.csv", index_col="id", dtype={"id": str})
    clin = pd.read_excel(DATA / "hsi" / "clinical.xlsx", sheet_name="Clinical Data")
    clin.index = clin["Case ID"].astype(str)
    h = s[s.cohort == "HistologyHSI-BC"].drop(columns=["time", "event"], errors="ignore").join(clin, how="inner")
    h["npi"] = 0.2 * h.Tumor_diameter / 10 + h.Tumor_histologic_grade + h.N.map({0: 1, 1: 2, 2: 3, 3: 3})
    for label, hh in [("HistologyHSI-BC HR+/HER2- (H&E)", h[((h.ER == 1) | (h.PR == 1)) & (h.HER2 == 0)]), ("HistologyHSI-BC all 47 (H&E)", h)]:
        x = pd.DataFrame({"score": std(hh.score), "npi": hh.npi})
        res[label] = evaluate(hh.DFS, hh.Relapse.astype(int),
                              {"npi": hh.npi, "score": hh.score, "score+npi": mods["score+npi"].predict_log_partial_hazard(x)}, "npi", rng)

    # TCGA HR+/HER2− (H&E 교차검증 점수, 등급 없음 → T · 림프절 · 나이)
    tc = s[s.cohort == "TCGA HR+/HER2- (CV)"].join(tcga_clinical())
    tc[CLIN] = tc[CLIN].fillna(tc[CLIN].median())
    x = tc[CLIN].assign(score=std(tc.score))
    res["TCGA HR+/HER2- (H&E CV)"] = evaluate(
        tc.time, tc.event.astype(int),
        {"clinical": mods["clinical"].predict_log_partial_hazard(tc[CLIN]), "score": tc.score,
         "score+clinical": mods["score+clinical"].predict_log_partial_hazard(x)}, "clinical", rng)
    OUT.write_text(json.dumps(res, indent=2, ensure_ascii=False))
    print(json.dumps(res, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
