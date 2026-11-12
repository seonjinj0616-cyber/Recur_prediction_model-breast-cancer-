"""C2 + 임상 결합 (방법 2): METABRIC에서 C2 점수와 임상 변수로 Cox를 학습하고, TCGA에서 H&E 점수를 넣어 평가한다.

임상 변수는 TCGA에도 있는 것만: T 병기(종양 크기로 T1 ≤ 20 mm / T2 ≤ 50 / T3 > 50, TCGA는 ajcc_tumor_pathologic_pt),
log(1 + 양성 림프절 수), 나이. 조직학적 등급은 TCGA에 없어 뺌.
C2 점수는 데이터마다 표준편차 1로 맞춤 (METABRIC: 측정 RNA 점수, TCGA: H&E 점수 gene · direct 또는 측정 RNA 점수).
학습: METABRIC HR+/HER2− 재발 기록 완전, 10년 원격 재발. 평가: TCGA 10년 PFI C-index (임상만 / C2만 / C2 + 임상).
먼저 c2_signature(METABRIC 점수)와 서버의 c2_he_direct(c2_he_scores.csv)를 실행하고 그 csv를 outputs/module1/에 둔다.
실행: python -m module1_recurrence.models.scripts.c2_clinical
결과: outputs/module1/c2_clinical.json
"""
import json

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

from module1_recurrence.data.scripts.external_fetch import DATA

OUT_DIR = DATA.parent / "outputs" / "module1"
TEN_YEARS = 3652.5
N_BOOT = 500
CLIN = ["t_stage", "log_nodes", "age"]


def t_from_size(mm):
    return np.select([mm <= 20, mm <= 50, mm > 50], [1, 2, 3], np.nan)


def tcga_clinical():
    c = pd.read_csv(DATA / "tcga" / "raw" / "nationwidechildrens.org_clinical_patient_brca.txt", sep="\t", skiprows=[1, 2])
    c = c.set_index("bcr_patient_barcode")
    t = c.ajcc_tumor_pathologic_pt.astype(str).str.extract(r"T(\d)")[0].astype(float).clip(upper=3)
    m = pd.read_csv(DATA / "tcga" / "manifest.csv").set_index("patient_id")
    return pd.DataFrame({"t_stage": t, "age": pd.to_numeric(c.age_at_diagnosis, errors="coerce")}).reindex(m.index).assign(
        log_nodes=np.log1p(pd.to_numeric(m.nodes, errors="coerce")))


def std(s):
    return (s - s.mean()) / s.std()


def main():
    rng = np.random.default_rng(0)
    mb = pd.read_csv(DATA / "metabric" / "cohort.csv")
    mb = mb[mb.complete_rec].set_index("patient_id")
    mb = mb.join(pd.read_csv(DATA / "metabric" / "c2_score.csv", index_col="patient_id").c2_raw)
    d = pd.DataFrame({"c2": std(mb.c2_raw), "t_stage": t_from_size(mb.tumor_size_mm), "log_nodes": np.log1p(mb.nodes),
                      "age": mb.age, "T": mb.dr_days.clip(upper=TEN_YEARS),
                      "E": ((mb.dr_event == 1) & (mb.dr_days <= TEN_YEARS)).astype(int)})
    d = d.fillna(d.median())
    models = {"clinical": CoxPHFitter().fit(d[[*CLIN, "T", "E"]], "T", "E"),
              "c2+clinical": CoxPHFitter().fit(d[["c2", *CLIN, "T", "E"]], "T", "E")}
    summary = {"metabric_hr": {k: m.hazard_ratios_.round(3).to_dict() for k, m in models.items()}}

    s = pd.read_csv(OUT_DIR / "c2_he_scores.csv", index_col=0)
    tc = s.join(tcga_clinical())
    tc[CLIN] = tc[CLIN].fillna(tc[CLIN].median())
    t, e = tc["T"].to_numpy(), tc["E"].to_numpy().astype(int)
    idx = [rng.integers(0, len(tc), len(tc)) for _ in range(N_BOOT)]
    scores = {"clinical only": models["clinical"].predict_log_partial_hazard(tc[CLIN]).to_numpy()}
    for k in ["truth", "gene", "direct"]:
        x = tc[CLIN].assign(c2=std(tc[k]))
        scores[f"{k} only"] = tc[k].to_numpy()
        scores[f"{k} + clinical"] = models["c2+clinical"].predict_log_partial_hazard(x).to_numpy()
    res = {}
    for k, sc in scores.items():
        boot = [concordance_index(t[i], -sc[i], e[i]) for i in idx]
        res[k] = {"c_index": round(concordance_index(t, -sc, e), 3), "95ci": [round(float(v), 3) for v in np.percentile(boot, [2.5, 97.5])]}
    summary["tcga_pfi"] = {"n": len(tc), "events_10y": int(e.sum()), **res}
    (OUT_DIR / "c2_clinical.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
