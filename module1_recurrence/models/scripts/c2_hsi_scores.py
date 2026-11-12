"""H&E → 원격 재발 평가용 점수 (서버): H&E → C2 재발 점수 모델을 만들고, 그 모델을 학습하지 않은 환자에게만 점수를 매긴다.

모듈 질문: HR+/HER2− 조기 유방암에서 H&E로 Oncotype처럼 10년 **원격** 재발 위험을 예측할 수 있는가 (평가는 원격 재발만).
입력은 CellViT 없이 전체 조직(썸네일 타일, thumb_tiles)만 쓴다 (2026-11-11 결정: 최종 점수 모델 · 평가는 CellViT 없이).
정답 점수 = TCGA 측정 RNA(log2 TPM, 학습 환자 안 z-score) × C2 가중치(METABRIC 학습, c2_signature_weights.csv) 합.
모델 (방법 1): 슬라이드 벡터(조직 타일 특징 평균, DX1) → 정답 점수 Ridge.
평가 대상
- TCGA HR+/HER2− (학습 코호트): 기관(TSS) 단위 5-fold 교차검증 점수 → 각 환자는 자기 병원을 학습에서 본 적 없는 모델의 점수.
  재발 = 10년 원격 전이 ({train}/distant.csv, 임상 new tumor event, manifest.load_distant). 측정 RNA 점수도 상한선으로 남김.
- HistologyHSI-BC: TCGA 전체로 학습한 모델의 점수 (재발 정보는 c2_hsi_eval에서 임상 표와 붙임).
실행: python -m module1_recurrence.models.scripts.c2_hsi_scores --train data/sj/tcga_thumb --hsi data/sj/hsi --weights data/sj/tcga/c2_he/c2_signature_weights.csv
결과: {hsi}/c2_validation_scores.csv (cohort, id, score, TCGA는 measured · recommended · time · event)
"""
import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from module1_recurrence.models.scripts.p_genes import ALPHAS

TEN_YEARS = 3652.5


def slide_means(feat_dir, pattern="*.h5"):
    vec = {}
    for f in sorted(feat_dir.glob(pattern)):
        with h5py.File(f, "r") as h:
            vec[f.stem] = h["features"][:].astype(np.float32).mean(0)
    return pd.DataFrame(vec).T


def model():
    return make_pipeline(StandardScaler(), RidgeCV(alphas=ALPHAS))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", type=Path, required=True, help="TCGA HR+/HER2− 폴더 (features/ · expr.parquet · manifest.csv · distant.csv)")
    ap.add_argument("--hsi", type=Path, required=True)
    ap.add_argument("--weights", type=Path, required=True)
    args = ap.parse_args()

    w = pd.read_csv(args.weights).set_index("covariate").weight
    x = slide_means(args.train / "features", "TCGA-*DX1.h5")
    x.index = x.index.str[:12]
    expr = pd.read_parquet(args.train / "expr.parquet")
    m = pd.read_csv(args.train / "manifest.csv").set_index("patient_id")
    dist = pd.read_csv(args.train / "distant.csv", index_col="patient_id")
    pids = x.index.intersection(expr.index).intersection(m.index)
    genes = [g for g in w.index if g in expr.columns]
    e = expr.loc[pids, genes]
    truth = ((e - e.mean()) / e.std()).fillna(0) @ w[genes]

    oof = cross_val_predict(model(), x.loc[pids].to_numpy(), truth.to_numpy(), cv=GroupKFold(5), groups=pids.str[5:7].to_numpy())
    d = dist.reindex(pids)
    tcga = pd.DataFrame({"cohort": "TCGA HR+/HER2- (CV)", "score": oof, "measured": truth.to_numpy(),
                         "recommended": m.loc[pids, "recommended"].to_numpy(),
                         "time": d.dm_days.clip(upper=TEN_YEARS).to_numpy(),
                         "event": ((d.dm_event == 1) & (d.dm_days <= TEN_YEARS)).astype(int).to_numpy()}, index=pids)
    print("TCGA", len(pids), "patients, CV r with RNA score", round(float(np.corrcoef(oof, truth)[0, 1]), 3),
          "| distant 10y", int(tcga.event.sum()))

    full = model().fit(x.loc[pids].to_numpy(), truth.to_numpy())
    h = slide_means(args.hsi / "features")
    hsi = pd.DataFrame({"cohort": "HistologyHSI-BC", "score": full.predict(h.to_numpy())}, index=h.index)
    out = pd.concat([tcga, hsi])
    out.index.name = "id"
    out.to_csv(args.hsi / "c2_validation_scores.csv")
    print("HistologyHSI-BC", len(hsi), "slides | alpha", full[-1].alpha_)


if __name__ == "__main__":
    main()
