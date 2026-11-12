"""C2 H&E 개선 시험: H&E 특징에서 C2 재발 점수를 얻는 방법들을 TCGA에서 같은 기준으로 비교한다.

정답 점수 = TCGA 측정 RNA(log2 TPM, 환자 안 z-score)에 C2 가중치(METABRIC 학습, c2_signature_weights.csv)를 곱한 합.
방법
- gene: 유전자 243개를 하나씩 예측(p_genes --save_pred, 기관 단위 교차검증)한 뒤 z-score → 같은 가중치로 합 (지금 방식)
- direct: 환자 벡터(맞춤 타일 평균, 1,536차원) → 점수 하나를 Ridge로 바로 예측 (기관 단위 5-fold, 교차검증 예측)
평가 (같은 환자): ① 정답 점수와의 Pearson r (점수 재현도, 주 지표) ② 10년 PFI C-index (사건이 적어 참고)
실행: python -m module1_recurrence.models.scripts.c2_he_direct --root data/sj/tcga
결과: {root}/c2_he/c2_he_direct.json, c2_he_scores.csv
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines.utils import concordance_index

from module1_recurrence.models.scripts.p_genes import ALPHAS, patient_vectors
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

TEN_YEARS = 3652.5
N_BOOT = 500


def z(df):
    return ((df - df.mean()) / df.std()).replace([np.inf, -np.inf], np.nan).fillna(0)


def direct(x, target, groups):
    model = make_pipeline(StandardScaler(), RidgeCV(alphas=ALPHAS))
    return cross_val_predict(model, x, target, cv=GroupKFold(5), groups=groups)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True, help="features/ · expr.parquet · manifest.csv · p_genes/ 가 있는 폴더")
    ap.add_argument("--weights", type=Path, default=None, help="c2_signature_weights.csv (기본: {root}/c2_he/)")
    args = ap.parse_args()
    out = args.root / "c2_he"
    out.mkdir(exist_ok=True)
    w = pd.read_csv(args.weights or out / "c2_signature_weights.csv").set_index("covariate").weight

    x, _, _ = patient_vectors(args.root, "matched", False, 20)
    pred = pd.read_parquet(args.root / "p_genes" / "p_genes_matched_dx1_tss_pred.parquet")
    measured = pd.read_parquet(args.root / "expr.parquet")
    m = pd.read_csv(args.root / "manifest.csv").set_index("patient_id")
    pids = x.index.intersection(pred.index).intersection(measured.index)
    genes = [g for g in w.index if g in pred.columns]
    w = w[genes]

    truth = z(measured.loc[pids, genes]) @ w
    s = pd.DataFrame({"truth": truth,
                      "gene": z(pred.loc[pids, genes]) @ w,
                      "direct": direct(x.loc[pids].to_numpy(), truth.to_numpy(), pids.str[5:7].to_numpy())}, index=pids)
    t = m.loc[pids, "PFI_time"].clip(upper=TEN_YEARS).to_numpy()
    e = ((m.loc[pids, "PFI"] == 1) & (m.loc[pids, "PFI_time"] <= TEN_YEARS)).to_numpy()
    s["T"], s["E"] = t, e
    s.to_csv(out / "c2_he_scores.csv")

    rng = np.random.default_rng(0)
    idx = [rng.integers(0, len(s), len(s)) for _ in range(N_BOOT)]
    res = {"n": len(s), "events_10y": int(e.sum()), "genes": len(genes)}
    for k in ["truth", "gene", "direct"]:
        ci = lambda i: concordance_index(t[i], -s[k].to_numpy()[i], e[i])  # noqa: E731
        boot = [ci(i) for i in idx]
        res[k] = {"r_with_truth": round(float(np.corrcoef(s[k], s.truth)[0, 1]), 3),
                  "pfi_c_index": round(ci(np.arange(len(s))), 3),
                  "pfi_95ci": [round(float(v), 3) for v in np.percentile(boot, [2.5, 97.5])]}
    (out / "c2_he_direct.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
