"""tile_sets.select: 고정 기준, RNA 조각 구성 맞춤(누적 종양 핵 비율이 목표에서 멈춤)."""
import numpy as np

from module1_recurrence.features.scripts.tile_sets import select

# 종양 핵 비율, 전체 핵 → 종양 핵 90, 80, 60, 40, 9, 0
RATIO = np.array([0.9, 0.8, 0.6, 0.4, 0.9, 0.0])
CELLS = np.array([100, 100, 100, 100, 10, 50])


def test_fixed():
    assert select(RATIO, CELLS, "tumor").tolist() == [True, True, False, False, False, False]
    assert select(RATIO, CELLS, "all").all()


def test_matched_stops_at_target():
    # 후보는 종양 핵 ≥ 50인 앞 3개. 누적 비율: 0.9 → 0.85 → 0.767
    assert select(RATIO, CELLS, "matched", purity=0.8).tolist() == [True, True, False, False, False, False]
    assert select(RATIO, CELLS, "matched", purity=0.7).tolist() == [True, True, True, False, False, False]
    assert select(RATIO, CELLS, "matched", purity=0.5).sum() == 3  # 후보를 다 담아도 목표보다 높으면 후보 전부


def test_matched_default_and_minimum():
    assert select(RATIO, CELLS, "matched").tolist() == select(RATIO, CELLS, "matched", purity=0.8).tolist()
    assert select(RATIO, CELLS, "matched", purity=0.99).sum() == 1  # 목표가 너무 높아도 가장 좋은 타일 1개
    assert select([0.1], [20], "matched", purity=0.8).sum() == 0  # 후보 없음
