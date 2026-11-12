"""외부 검증용 H&E 점수 (서버): TCGA HR+/HER2− 513명으로 H&E → C2 점수 모델을 학습하고, 학습에 쓰지 않은 환자에만 점수를 매긴다.

CellViT 없이 전체 조직(썸네일 타일, thumb_tiles)만 쓴다 (2026-11-11 결정: 최종 점수 모델 · 검증은 CellViT 없이).
정답 점수 = TCGA 측정 RNA(log2 TPM, 학습 환자 안 z-score) × C2 가중치(METABRIC 학습, c2_signature_weights.csv) 합.
모델 (방법 1): 슬라이드 벡터(조직 타일 특징 평균, DX1) → 정답 점수 Ridge, 학습 환자 전체로 한 번 학습.
검증 (학습에 쓰지 않음):
- HistologyHSI-BC 47명
- TCGA 다른 아형 (data/tcga/other, 학습 · 유전자 선택 어디에도 쓰지 않음), 진단 시 원격 전이(M1) 제외.
  재발 = 10년 원격 전이 (임상 new tumor event, manifest.load_distant). PFI(국소 재발 · 진행 포함)도 참고로 남김.
  측정 RNA 점수(같은 학습 평균 · 표준편차로 z-score)도 같이 남겨 상한선으로 비교
실행: python -m module1_recurrence.models.scripts.c2_hsi_scores --train data/sj/tcga_thumb --hsi data/sj/hsi --other data/sj/tcga_other --weights data/sj/tcga/c2_he/c2_signature_weights.csv
결과: {hsi}/c2_validation_scores.csv (cohort, id, score, TCGA 다른 아형은 measured · group · 원격 전이 시간 · 사건 · PFI)
"""
import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
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


def tcga_vectors(root):
    x = slide_means(root / "features", "TCGA-*DX1.h5")
    x.index = x.index.str[:12]
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", type=Path, required=True, help="학습 TCGA 폴더 (features/ · expr.parquet · manifest.csv)")
    ap.add_argument("--hsi", type=Path, required=True)
    ap.add_argument("--other", type=Path, required=True, help="TCGA 다른 아형 폴더 (features/ · expr.parquet · manifest.csv)")
    ap.add_argument("--weights", type=Path, required=True)
    args = ap.parse_args()

    w = pd.read_csv(args.weights).set_index("covariate").weight
    x = tcga_vectors(args.train)
    expr = pd.read_parquet(args.train / "expr.parquet")
    pids = x.index.intersection(expr.index)
    genes = [g for g in w.index if g in expr.columns]
    mu, sd = expr.loc[pids, genes].mean(), expr.loc[pids, genes].std()
    rna_score = lambda e: (((e[genes] - mu) / sd).replace([np.inf, -np.inf], np.nan).fillna(0) @ w[genes])  # noqa: E731
    model = make_pipeline(StandardScaler(), RidgeCV(alphas=ALPHAS)).fit(x.loc[pids].to_numpy(), rna_score(expr.loc[pids]).to_numpy())
    print("trained on", len(pids), "TCGA HR+/HER2- patients | alpha", model[-1].alpha_)

    h = slide_means(args.hsi / "features")
    out = [pd.DataFrame({"cohort": "HistologyHSI-BC", "score": model.predict(h.to_numpy())}, index=h.index)]

    xo = tcga_vectors(args.other)
    if len(xo):
        m = pd.read_csv(args.other / "manifest.csv").set_index("patient_id")
        m = m[~m.m1 & m.dm_days.notna()]  # 원격 재발 검증: 진단 시 원격 전이(M1) · 추적 시간 없는 환자 제외
        eo = pd.read_parquet(args.other / "expr.parquet")
        po = xo.index.intersection(m.index).intersection(eo.index)
        assert not set(po) & set(pids), "검증 환자가 학습에 들어 있음"
        out.append(pd.DataFrame({"cohort": "TCGA other subtypes", "score": model.predict(xo.loc[po].to_numpy()),
                                 "measured": rna_score(eo.loc[po]).to_numpy(), "group": m.loc[po, "group"].to_numpy(),
                                 "time": m.loc[po, "dm_days"].clip(upper=TEN_YEARS).to_numpy(),
                                 "event": ((m.loc[po, "dm_event"] == 1) & (m.loc[po, "dm_days"] <= TEN_YEARS)).astype(int).to_numpy(),
                                 "pfi_time": m.loc[po, "PFI_time"].clip(upper=TEN_YEARS).to_numpy(),
                                 "pfi_event": ((m.loc[po, "PFI"] == 1) & (m.loc[po, "PFI_time"] <= TEN_YEARS)).astype(int).to_numpy()},
                                index=po))
        print("TCGA other subtypes scored", len(po))
    res = pd.concat(out)
    res.index.name = "id"
    res.to_csv(args.hsi / "c2_validation_scores.csv")
    print(res.groupby("cohort").size())


if __name__ == "__main__":
    main()
