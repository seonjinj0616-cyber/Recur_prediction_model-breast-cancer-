"""A5: 비교 기준선 점수를 METABRIC · SCAN-B 환자마다 계산하고, 같은 환자에서 10년 원격 재발 C-index를 비교한다.

- pseudo_rs: A3 결과 (labels.scripts.pseudo_rs_prognosis)
- proliferation: Oncotype 증식 유전자 5개 z-score 평균 (HR+/HER2− 전체 안에서 z-score)
- npi: Nottingham Prognostic Index = 0.2 × 크기(cm) + 등급 + 림프절 단계(0개 1, 1–3개 2, 4개 이상 3).
  METABRIC은 제공 열, SCAN-B는 Size.mm · NHG · LN.spec으로 계산
- ror: PAM50 ROR (SCAN-B 제공 NCN.ROR.asT0. asT1은 크기 항만 상수로 달라 순위가 같다). METABRIC은 제공값이 없어 비움
- sequoia: SEQUOIA 272개 유전자 서명 (Pizurica et al. 2024, Supplementary Data 5 계수, TCGA log2(FPKM-UQ + 1)에서 LASSO Cox 학습).
  파일(CC BY-NC-ND 4.0, 저장소에 올리지 않음): https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-54182-5/MediaObjects/41467_2024_54182_MOESM4_ESM.xlsx
  → data/sequoia/raw/sequoia_supplementary_data_1-5.xlsx
  점수 = Σ 계수 × (log2 발현 − 코호트 평균). 공통 유전자에 없는 13개는 빼고 259개로 계산
  (계수가 표준화하지 않은 log 값 기준이라 z-score 대신 평균만 뺀다. 마이크로어레이에서는 플랫폼 차이로 근사)
CTS5와 무작위 유전자 세트는 여기서 계산하지 않는다 (CTS5 공식 원문 확인 전, 무작위 세트는 C2에서 최종 서명과 같은 개수로).
실행: python -m module1_recurrence.scripts.baselines
결과: data/{metabric,scanb}/baselines.csv, outputs/module1/a5_baselines.json
"""
import json

import numpy as np
import pandas as pd
from lifelines.utils import concordance_index

from module1_recurrence.data.scripts.expression import load, zscore
from module1_recurrence.data.scripts.external_fetch import DATA, SCANB
from module1_recurrence.labels.scripts.pseudo_rs import PROLIFERATION
from module1_recurrence.scripts.cohort import RECOMMENDED

OUT_DIR = DATA.parent / "outputs" / "module1"
TEN_YEARS = 3652.5
N_BOOT = 200
SEQUOIA_FILE = DATA / "sequoia" / "raw" / "sequoia_supplementary_data_1-5.xlsx"
NPI_NODE_STAGE = {"N0": 1, "SubMicroMet": 1, "1to3": 2, "4toX": 3}


def sequoia_coef():
    s = pd.read_excel(SEQUOIA_FILE, sheet_name="Supp. Data 5", header=2)
    return s.set_index("Gene")["Hazard ratio (Coeff.)"]


def scores(name, cohort):
    log_expr = load(name)
    expr = zscore(log_expr)
    out = pd.DataFrame(index=cohort.patient_id)
    out["pseudo_rs"] = pd.read_csv(DATA / name / "pseudo_rs.csv", index_col="patient_id").pseudo_rs
    out["proliferation"] = expr[PROLIFERATION].mean(axis=1)
    coef = sequoia_coef()
    coef = coef[coef.index.isin(log_expr.columns)]
    centered = log_expr[coef.index] - log_expr[coef.index].mean()
    out["sequoia"] = centered.fillna(0) @ coef
    c = cohort.set_index("patient_id")
    if name == "metabric":
        out["npi"] = c.npi
    else:
        clin = pd.read_excel(SCANB / "scanb_clinical_2023-01-13.xlsx", sheet_name="SCANB.9206")
        clin = clin[clin["Follow.up.cohort"]].set_index("Patient").reindex(c.index)
        out["npi"] = 0.2 * clin["Size.mm"] / 10 + clin.NHG + clin["LN.spec"].map(NPI_NODE_STAGE)
        out["ror"] = clin["NCN.ROR.asT0"]
    return out


def compare(df, cols, rng):
    """cols 점수가 모두 있는 같은 환자에서 C-index와 부트스트랩 95% CI."""
    df = df.dropna(subset=cols)
    t = df.dr_days.clip(upper=TEN_YEARS)
    e = ((df.dr_event == 1) & (df.dr_days <= TEN_YEARS)).astype(int)
    idx = [rng.integers(0, len(df), len(df)) for _ in range(N_BOOT)]
    res = {"n": len(df), "events_10y": int(e.sum())}
    for c in cols:
        s = -df[c].to_numpy()  # 점수가 높을수록 빨리 재발
        boot = [concordance_index(t.iloc[i], s[i], e.iloc[i]) for i in idx]
        res[c] = {"c_index": round(concordance_index(t, s, e), 3),
                  "c_index_95ci": [round(float(x), 3) for x in np.percentile(boot, [2.5, 97.5])]}
    return res


def main():
    rng = np.random.default_rng(0)
    results = {}
    for name in ["metabric", "scanb"]:
        cohort = pd.read_csv(DATA / name / "cohort.csv")
        s = scores(name, cohort)
        s.round(4).to_csv(DATA / name / "baselines.csv")

        df = cohort.drop(columns="npi", errors="ignore").dropna(subset=["dr_event", "dr_days"]).join(s, on="patient_id")
        if name == "metabric":
            df = df[df.complete_rec]
        cols = list(s.columns[s.notna().any()])
        results[name] = {"all HR+/HER2-": compare(df, cols, rng),
                         "recommended": compare(df[df.group.isin(RECOMMENDED)], cols, rng)}

        print(f"\n{name}")
        for label, r in results[name].items():
            print(f"  {label}: n={r['n']}, 10년 원격 재발 {r['events_10y']}건")
            for c in cols:
                print(f"    {c:14s} C-index {r[c]['c_index']} {r[c]['c_index_95ci']}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "a5_baselines.json").write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print("\nsaved", OUT_DIR / "a5_baselines.json")


if __name__ == "__main__":
    main()
