"""C2 + 임상 결합: 재발 점수에 임상 정보를 더하면 원격 재발 예측이 얼마나 좋아지는지 본다.

결합 가중치는 METABRIC HR+/HER2− 1,251명(10년 원격 재발)에서 Cox로 학습해 고정하고, 평가 데이터에서는 다시 맞추지 않는다.
점수는 데이터마다 표준편차 1로 맞춘다 (METABRIC: 측정 RNA C2, 평가: RNA 또는 H&E 점수).
- NPI 결합 (SCAN-B, HistologyHSI-BC): Cox(점수 + NPI). NPI = 0.2 × 크기(cm) + 등급 + 림프절 단계.
- 임상 결합 (TCGA, 등급이 없어 NPI 대신): Cox(점수 + T 병기 + log(1 + 양성 림프절) + 나이) — c2_clinical과 같은 변수.
- CTS5 결합 (모든 데이터): Cox(점수 + CTS5). CTS5 = 0.438 × 림프절 + 0.988 × (0.093 × 크기 − 0.001 × 크기² + 0.375 × 등급 + 0.017 × 나이)
  (Dowsett et al. JCO 2018, 원문 확인: 크기 mm · 30 mm에서 자름, 림프절 0 / 1 / 2–3 / 4–9 / 10개 이상 → 0–4).
  데이터별 입력과 근사
  · METABRIC: 양성 림프절 수 · 크기 mm · 등급 · 나이 그대로
  · SCAN-B: 림프절이 범주(N0 · SubMicroMet → 0, 1to3 → 1.5(1과 2–3의 가운데), 4toX → 3), 나이는 5세 단위(age_range)
  · HistologyHSI-BC: 양성 림프절 = LN_MICRO_number + LN_MACRO_number (ITC 제외), 지름 mm
  · TCGA: 크기가 pT만 있어 T1a 3 · T1b 7.5 · T1/T1c 15 · T2 이상 30 mm(자르는 값)로 근사. 림프절 수가 없으면 pN으로
    (N0 → 0, N1 → 1.5, N2 → 3, N3 → 4, NX → 중앙값). 등급은 data/tcga/grade.csv (판독 + 직접 지정 1명 + Y-05 14명).
    주 결과는 513명, 민감도 분석은 판독 등급(thennavan2021)이 있는 환자만.
평가: 임상만 · 점수만 · 점수 + 임상의 10년 원격 재발 C-index와, 같은 부트스트랩 표본에서 (점수 + 임상) − (임상만) 차이의 95% CI.
- SCAN-B 검사 권고 (측정 RNA C2)
- HistologyHSI-BC HR+/HER2− (H&E 점수, 주 평가) · 47명 전체(보조)
- TCGA HR+/HER2− (H&E 교차검증 점수)
먼저 data/scripts/tcga_grade(data/tcga/grade.csv), c2_signature(c2_score.csv), c2_hsi_scores · c2_hsi_eval(data/hsi/c2_validation_scores.csv)을 실행한다.
실행: python -m module1_recurrence.models.scripts.c2_combined
결과: outputs/module1/c2_combined.json
"""
import json

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

from module1_recurrence.data.scripts.external_fetch import DATA
from module1_recurrence.data.scripts.gdc_fetch import RAW
from module1_recurrence.models.scripts.c2_clinical import CLIN, std, t_from_size, tcga_clinical
from module1_recurrence.scripts.cohort import RECOMMENDED

OUT = DATA.parent / "outputs" / "module1" / "c2_combined.json"
TEN_YEARS = 3652.5
N_BOOT = 1000
CTS5_PT_MM = {"T1a": 3, "T1b": 7.5, "T1": 15, "T1c": 15}  # 그 밖(T2 이상)은 30 mm
CTS5_PN = {"N0": 0, "N1": 1.5, "N2": 3, "N3": 4}


def node_cat(n):
    """양성 림프절 수 → CTS5 범주 0 / 1 / 2–3 / 4–9 / 10개 이상 = 0–4."""
    return pd.cut(n, [-np.inf, 0, 1, 3, 9, np.inf], labels=False)


def cts5(size_mm, nodes, grade, age):
    """nodes는 CTS5 범주(0–4) 값."""
    s = np.minimum(size_mm, 30)
    return 0.438 * nodes + 0.988 * (0.093 * s - 0.001 * s ** 2 + 0.375 * grade + 0.017 * age)


def tcga_cts5():
    """TCGA HR+/HER2− CTS5와 등급 출처."""
    c = pd.read_csv(RAW / "nationwidechildrens.org_clinical_patient_brca.txt", sep="\t", skiprows=[1, 2]).set_index("bcr_patient_barcode")
    m = pd.read_csv(DATA / "tcga" / "manifest.csv").set_index("patient_id")
    g = pd.read_csv(DATA / "tcga" / "grade.csv", index_col="patient_id")
    c = c.reindex(m.index)
    pt = c.ajcc_tumor_pathologic_pt.astype(str)
    size = pt.map(CTS5_PT_MM).where(pt.str.startswith("T1"), 30.0)
    pn = c.ajcc_nodes_pathologic_pn.astype(str).str.extract(r"^(N\d)")[0].map(CTS5_PN)
    nodes = node_cat(pd.to_numeric(m.nodes, errors="coerce")).astype(float).fillna(pn)
    nodes = nodes.fillna(nodes.median())
    age = pd.to_numeric(c.age_at_diagnosis, errors="coerce")
    return pd.DataFrame({"cts5": cts5(size, nodes, g.grade.reindex(m.index), age.fillna(age.median())),
                         "grade_source": g.source.reindex(m.index)})


def metabric_models():
    mb = pd.read_csv(DATA / "metabric" / "cohort.csv")
    mb = mb[mb.complete_rec].set_index("patient_id").join(pd.read_csv(DATA / "metabric" / "c2_score.csv", index_col="patient_id").c2_raw)
    d = pd.DataFrame({"score": std(mb.c2_raw), "npi": mb.npi, "t_stage": t_from_size(mb.tumor_size_mm), "log_nodes": np.log1p(mb.nodes),
                      "age": mb.age, "cts5": cts5(mb.tumor_size_mm, node_cat(mb.nodes), mb.grade, mb.age), "T": mb.dr_days.clip(upper=TEN_YEARS),
                      "E": ((mb.dr_event == 1) & (mb.dr_days <= TEN_YEARS)).astype(int)})
    d = d.fillna(d.median())
    fit = lambda cols: CoxPHFitter().fit(d[[*cols, "T", "E"]], "T", "E")  # noqa: E731
    return {"npi": fit(["npi"]), "score+npi": fit(["score", "npi"]), "clinical": fit(CLIN), "score+clinical": fit(["score", *CLIN]),
            "cts5": fit(["cts5"]), "score+cts5": fit(["score", "cts5"])}


def evaluate(t, e, scores, rng):
    """C-index와, 같은 부트스트랩 표본에서 (score+X) − X 차이."""
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
            base = k.split("+", 1)[1]
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
        pd.read_csv(DATA / "scanb" / "baselines.csv", index_col="patient_id").npi)
    d["cts5"] = cts5(d.tumor_size_mm, d.nodes.map({"N0": 0, "SubMicroMet": 0, "1to3": 1.5, "4toX": 3}), d.grade, d.age_range)
    d = d.dropna(subset=["c2_raw", "npi", "cts5"])
    x = pd.DataFrame({"score": std(d.c2_raw), "npi": d.npi, "cts5": d.cts5})
    res["SCAN-B recommended (RNA)"] = evaluate(
        d.dr_days.clip(upper=TEN_YEARS), ((d.dr_event == 1) & (d.dr_days <= TEN_YEARS)).astype(int),
        {"npi": d.npi, "cts5": d.cts5, "score": d.c2_raw, "score+npi": mods["score+npi"].predict_log_partial_hazard(x),
         "score+cts5": mods["score+cts5"].predict_log_partial_hazard(x)}, rng)

    # HistologyHSI-BC (H&E 점수)
    s = pd.read_csv(DATA / "hsi" / "c2_validation_scores.csv", index_col="id", dtype={"id": str})
    clin = pd.read_excel(DATA / "hsi" / "clinical.xlsx", sheet_name="Clinical Data")
    clin.index = clin["Case ID"].astype(str)
    h = s[s.cohort == "HistologyHSI-BC"].drop(columns=["time", "event"], errors="ignore").join(clin, how="inner")
    h["npi"] = 0.2 * h.Tumor_diameter / 10 + h.Tumor_histologic_grade + h.N.map({0: 1, 1: 2, 2: 3, 3: 3})
    h["cts5"] = cts5(h.Tumor_diameter, node_cat(h.LN_MICRO_number + h.LN_MACRO_number), h.Tumor_histologic_grade,
                     h["Age at Diagnosis (Years)"])
    for label, hh in [("HistologyHSI-BC HR+/HER2- (H&E)", h[((h.ER == 1) | (h.PR == 1)) & (h.HER2 == 0)]), ("HistologyHSI-BC all 47 (H&E)", h)]:
        x = pd.DataFrame({"score": std(hh.score), "npi": hh.npi, "cts5": hh.cts5})
        res[label] = evaluate(hh.DFS, hh.Relapse.astype(int),
                              {"npi": hh.npi, "cts5": hh.cts5, "score": hh.score,
                               "score+npi": mods["score+npi"].predict_log_partial_hazard(x),
                               "score+cts5": mods["score+cts5"].predict_log_partial_hazard(x)}, rng)

    # TCGA HR+/HER2− (H&E 교차검증 점수, 등급 없음 → T · 림프절 · 나이)
    tc = s[s.cohort == "TCGA HR+/HER2- (CV)"].join(tcga_clinical()).join(tcga_cts5())
    tc[CLIN] = tc[CLIN].fillna(tc[CLIN].median())
    for label, tt in [("TCGA HR+/HER2- (H&E CV)", tc), ("TCGA HR+/HER2- read grade only (H&E CV)", tc[tc.grade_source == "thennavan2021"])]:
        x = tt[CLIN].assign(score=std(tt.score), cts5=tt.cts5)
        res[label] = evaluate(
            tt.time, tt.event.astype(int),
            {"clinical": mods["clinical"].predict_log_partial_hazard(tt[CLIN]), "cts5": tt.cts5, "score": tt.score,
             "score+clinical": mods["score+clinical"].predict_log_partial_hazard(x),
             "score+cts5": mods["score+cts5"].predict_log_partial_hazard(x)}, rng)
    OUT.write_text(json.dumps(res, indent=2, ensure_ascii=False))
    print(json.dumps(res, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
