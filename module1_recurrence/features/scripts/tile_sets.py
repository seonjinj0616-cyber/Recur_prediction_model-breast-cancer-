"""B4에서 환자마다 쓸 타일을 고른다. 입력은 features h5의 tumor_ratio · n_cells (환자의 슬라이드를 이어 붙인 것).

h5의 tumor_tile 표시는 추출할 때의 기준(0.5)이라 쓰지 않고, 여기서 핵 수로 다시 판정한다.
종양 핵 수 = round(tumor_ratio × n_cells).

- all: 조직 타일 전부 (비교용)
- tumor: 고정 기준. 종양 핵 비율 ≥ 0.7 & 종양 핵 ≥ 50 (128 µm 타일당)
- matched: RNA 조각 구성 맞춤. 종양 핵 ≥ 50인 타일을 비율 높은 순으로 담다가,
  담은 타일 전체의 종양 핵 비율(종양 핵 합 / 전체 핵 합)이 목표값까지 내려오면 멈춘다.
  목표값은 학습 때 환자의 RNA 조각 종양 핵 % (data/tcga/rna_purity.csv), 예측 때 코호트 중앙값 80%.
  비율이 높은 타일부터 담으므로 종양 중심에서 경계 쪽으로 넓혀 가는 macrodissection과 비슷하다.
"""
import numpy as np

MIN_RATIO = 0.7
MIN_NEOPLASTIC = 50
DEFAULT_PURITY = 0.8  # TCGA HR+/HER2− 515명 RNA 조각 종양 핵 % 중앙값


def select(tumor_ratio, n_cells, mode, purity=None, min_neoplastic=MIN_NEOPLASTIC):
    """고를 타일의 bool 마스크. purity는 0–1 (없으면 DEFAULT_PURITY)."""
    tumor_ratio, n_cells = np.asarray(tumor_ratio, dtype=float), np.asarray(n_cells)
    neo = np.rint(tumor_ratio * n_cells)
    if mode == "all":
        return np.ones(len(n_cells), dtype=bool)
    if mode == "tumor":
        return (tumor_ratio >= MIN_RATIO) & (neo >= min_neoplastic)
    if mode != "matched":
        raise ValueError(mode)

    target = DEFAULT_PURITY if purity is None or np.isnan(purity) else purity
    cand = np.flatnonzero(neo >= min_neoplastic)
    order = cand[np.argsort(-tumor_ratio[cand], kind="stable")]
    # 비율 높은 순으로 더하므로 누적 비율은 줄어들기만 한다 → 목표 이상인 앞부분까지 담는다
    cum = np.cumsum(neo[order]) / np.cumsum(n_cells[order])
    k = max(1, int((cum >= target).sum())) if len(order) else 0
    keep = np.zeros(len(n_cells), dtype=bool)
    keep[order[:k]] = True
    return keep
