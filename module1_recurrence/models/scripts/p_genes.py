"""B4 (단계 2): H&E 특징으로 유전자 발현을 예측하고 SEQUOIA 기준으로 "잘 예측된 유전자"(P 목록)를 고른다.

1. 환자마다 슬라이드 h5(features · tumor_ratio · n_cells)를 이어 붙이고 tile_sets.select로 타일을 고른 뒤 평균 → 환자 벡터
   (--tiles matched: RNA 조각 종양 핵 %에 맞춤 / tumor: 고정 기준 / all: 조직 타일 전부)
   유전자는 코호트 전체(515명) 중 절반 이상에서 TPM ≥ 1인 것만 쓴다 (12,758 / 17,813개). TPM 0은 빠진 값이 아니라
   "검출 안 됨"이라 값을 빼지 않고, 대부분 환자에서 측정되지 않는 유전자를 통째로 뺀다.
2. 기관(TSS) 단위 5-fold (같은 병원 환자는 한 fold에, `--cv random`이면 환자 단위 무작위),
   StandardScaler + RidgeCV(유전자마다 alpha) → 유전자 log2(TPM + 1)의 교차검증 예측.
   병원마다 스캐너 · 염색 · RNA 처리 인력이 달라, 무작위로 나누면 병원 차이를 외운 것도 예측력처럼 보이기 때문.
3. 대조: 같은 차원의 무작위 특징으로 같은 과정 (SEQUOIA의 "학습하지 않은 무작위 모델"에 해당)
4. 유전자마다 SEQUOIA 기준 세 가지를 모두 만족하면 잘 예측됨:
   ① Pearson r > 0, P < 0.05  ② 대조보다 r이 유의하게 큼 (Steiger · Meng 검정, 한쪽 P < 0.05, BH < 0.2)
   ③ RMSE가 대조보다 작음.  순위는 Pearson r 내림차순.

먼저 data.scripts.tcga_expression으로 data/tcga/expr.parquet, rna_purity로 rna_purity.csv를 만들어 {root}에 둔다.
실행: python -m module1_recurrence.models.scripts.p_genes --root data/sj/tcga [--tiles matched] [--all_slides]
결과: {out}/p_genes_{tiles}_{dx1|all}_tss.csv (random이면 _tss 없음) (유전자별 r · P · 대조 r · Steiger · RMSE · well), 같은 이름 .json 요약
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupKFold, KFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from module1_recurrence.data.scripts.gdc_fetch import RAW
from module1_recurrence.features.scripts.tile_sets import select

ALPHAS = np.logspace(0, 6, 13)


def patient_vectors(root, mode, all_slides, min_tiles):
    """환자 ID → 고른 타일 특징 평균. 고른 타일이 min_tiles개 미만인 환자는 뺀다."""
    import h5py  # 서버(특징이 있는 곳)에서만 필요
    from module1_recurrence.features.scripts.extract import ENCODER
    purity = (pd.read_csv(root / "rna_purity.csv").set_index("patient_id").tumor_nuclei / 100).to_dict()
    files = sorted((root / "features").glob("TCGA-*.h5"))
    if not all_slides:
        files = [f for f in files if f.stem.endswith("DX1")]
    by_patient = {}
    for f in files:
        by_patient.setdefault(f.stem[:12], []).append(f)

    vec, n_tiles, dropped = {}, {}, []
    for pid, paths in by_patient.items():
        parts = []
        for p in paths:
            with h5py.File(p, "r") as h:
                if h.attrs.get("encoder") != ENCODER:
                    continue
                parts.append((h["features"][:], h["tumor_ratio"][:], h["n_cells"][:]))
        if not parts:
            continue
        x, ratio, cells = (np.concatenate(a) for a in zip(*parts))
        keep = select(ratio, cells, mode, purity=purity.get(pid))
        if keep.sum() < min_tiles:
            dropped.append((pid, int(keep.sum())))
            continue
        vec[pid] = x[keep].astype(np.float32).mean(0)
        n_tiles[pid] = int(keep.sum())
    return pd.DataFrame(vec).T, pd.Series(n_tiles), dropped


def expressed_genes(expr, min_frac):
    """코호트 전체에서 TPM ≥ 1 (log2(TPM + 1) ≥ 1)인 환자 비율이 min_frac 이상인 유전자."""
    return expr.columns[(expr >= 1).mean() >= min_frac]


def cv_predict(x, y, seed, groups=None):
    """groups가 있으면 같은 그룹(기관)이 학습 · 평가에 나뉘지 않게 fold를 나눈다."""
    model = make_pipeline(StandardScaler(), RidgeCV(alphas=ALPHAS, alpha_per_target=True))
    cv = KFold(5, shuffle=True, random_state=seed) if groups is None else GroupKFold(5)
    return cross_val_predict(model, x, y, cv=cv, groups=groups)


def pearson_cols(a, b):
    ac, bc = a - a.mean(0), b - b.mean(0)
    return (ac * bc).sum(0) / np.sqrt((ac ** 2).sum(0) * (bc ** 2).sum(0) + 1e-12)


def steiger(r_y1, r_y2, r_12, n):
    """같은 y에 대한 두 상관 r(y, 예측1)과 r(y, 예측2)의 차이 (Meng, Rosenthal & Rubin 1992). 한쪽 P (1 > 2)."""
    r_y1, r_y2 = np.clip(r_y1, -0.999999, 0.999999), np.clip(r_y2, -0.999999, 0.999999)
    rm2 = (r_y1 ** 2 + r_y2 ** 2) / 2
    f = np.minimum((1 - r_12) / (2 * (1 - rm2)), 1)
    h = (1 - f * rm2) / (1 - rm2)
    z = (np.arctanh(r_y1) - np.arctanh(r_y2)) * np.sqrt((n - 3) / (2 * (1 - r_12) * h))
    return z, stats.norm.sf(z)


def bh(p):
    p = np.asarray(p)
    order = np.argsort(p)
    q = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1)
    return out


def judge(y, pred, pred_null):
    """유전자별 SEQUOIA 기준 표 (y · pred · pred_null: 환자 × 유전자)."""
    n = len(y)
    r, r_null, r_pp = pearson_cols(y, pred), pearson_cols(y, pred_null), pearson_cols(pred, pred_null)
    t = r * np.sqrt((n - 2) / np.maximum(1 - r ** 2, 1e-12))
    p = 2 * stats.t.sf(np.abs(t), n - 2)
    z, p_st = steiger(r, r_null, r_pp, n)
    rmse = np.sqrt(((y - pred) ** 2).mean(0))
    rmse_null = np.sqrt(((y - pred_null) ** 2).mean(0))
    df = pd.DataFrame({"r": r, "p": p, "r_null": r_null, "steiger_z": z, "steiger_p": p_st,
                       "steiger_fdr": bh(p_st), "rmse": rmse, "rmse_null": rmse_null})
    df["well"] = (df.r > 0) & (df.p < 0.05) & (df.steiger_p < 0.05) & (df.steiger_fdr < 0.2) & (df.rmse < df.rmse_null)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=RAW.parent, help="features/ · rna_purity.csv · expr.parquet 폴더 (서버: data/sj/tcga)")
    ap.add_argument("--tiles", default="matched", choices=["matched", "tumor", "all"])
    ap.add_argument("--all_slides", action="store_true", help="환자의 DX 슬라이드 전부 (기본: DX1만)")
    ap.add_argument("--min_tiles", type=int, default=20, help="고른 타일이 이보다 적은 환자는 뺀다")
    ap.add_argument("--min_expr_frac", type=float, default=0.5, help="TPM ≥ 1인 환자가 이 비율 이상인 유전자만")
    ap.add_argument("--out", type=Path, default=None, help="결과 폴더 (기본: {root}/p_genes)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--cv", default="tss", choices=["tss", "random"],
                    help="tss(기본): 같은 기관(TSS, 환자 ID 6–7번째 글자) 환자를 한 fold에 묶음. random: 환자 단위 무작위 (SEQUOIA 방식, 비교용)")
    ap.add_argument("--save_pred", action="store_true", help="교차검증 예측 발현(환자 × 유전자)을 {name}_pred.parquet로 저장 (C2 H&E 평가용)")
    args = ap.parse_args()

    out = args.out or args.root / "p_genes"
    out.mkdir(parents=True, exist_ok=True)
    x, n_tiles, dropped = patient_vectors(args.root, args.tiles, args.all_slides, args.min_tiles)
    expr = pd.read_parquet(args.root / "expr.parquet")
    expr = expr[expressed_genes(expr, args.min_expr_frac)]
    pids = x.index.intersection(expr.index)
    x, y = x.loc[pids].to_numpy(), expr.loc[pids]
    y = y.loc[:, y.std() > 0]
    print(f"patients {len(pids)} (dropped {len(dropped)} with < {args.min_tiles} tiles), genes {y.shape[1]}", flush=True)

    groups = pids.str[5:7].to_numpy() if args.cv == "tss" else None
    pred = cv_predict(x, y.to_numpy(), args.seed, groups)
    null = np.random.default_rng(args.seed).normal(size=x.shape)
    pred_null = cv_predict(null, y.to_numpy(), args.seed, groups)
    res = judge(y.to_numpy(), pred, pred_null)
    res.insert(0, "gene", y.columns)
    res = res.sort_values("r", ascending=False).reset_index(drop=True)
    res.insert(1, "rank", np.arange(1, len(res) + 1))

    name = f"p_genes_{args.tiles}_{'all' if args.all_slides else 'dx1'}" + ("_tss" if args.cv == "tss" else "")
    res.to_csv(out / f"{name}.csv", index=False)
    if args.save_pred:
        pd.DataFrame(pred.astype(np.float32), index=y.index, columns=y.columns).to_parquet(out / f"{name}_pred.parquet")
    summary = {"tiles": args.tiles, "all_slides": args.all_slides, "cv": args.cv, "patients": len(pids),
               "dropped": dropped, "median_tiles": int(n_tiles.loc[pids].median()), "min_expr_frac": args.min_expr_frac, "genes": len(res),
               "well_predicted": int(res.well.sum()), "median_r": round(float(res.r.median()), 3),
               "median_r_null": round(float(res.r_null.median()), 3),
               "well_r_ge_0.3": int((res.well & (res.r >= 0.3)).sum())}
    (out / f"{name}.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
