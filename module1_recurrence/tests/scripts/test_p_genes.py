"""p_genes: Steiger 검정, BH, SEQUOIA 기준 판정."""
import numpy as np

import pandas as pd

from module1_recurrence.models.scripts.p_genes import bh, expressed_genes, judge, steiger


def test_steiger_direction_and_equal():
    z, p = steiger(np.array([0.6]), np.array([0.1]), np.array([0.1]), 100)
    assert z[0] > 0 and p[0] < 0.001
    z, p = steiger(np.array([0.3]), np.array([0.3]), np.array([0.5]), 100)
    assert abs(z[0]) < 1e-9 and abs(p[0] - 0.5) < 1e-9


def test_bh_matches_definition():
    q = bh([0.01, 0.04, 0.03, 0.5])
    # 정렬 p 0.01, 0.03, 0.04, 0.5 → ×4/순위 0.04, 0.06, 0.0533, 0.5 → 뒤에서부터 최솟값 0.04, 0.0533, 0.0533, 0.5
    assert np.allclose(q, [0.04, 0.16 / 3, 0.16 / 3, 0.5])


def test_judge_separates_signal_from_noise():
    rng = np.random.default_rng(0)
    n = 200
    y = rng.normal(size=(n, 2))
    pred = np.c_[y[:, 0] + rng.normal(scale=0.5, size=n), rng.normal(size=n)]  # 0번만 신호
    null = rng.normal(size=(n, 2))
    res = judge(y, pred, null)
    assert res.well.tolist() == [True, False]


def test_expressed_genes_keeps_majority_expressed():
    # log2(TPM + 1): 1 이상이면 TPM ≥ 1. a는 4명 중 3명, b는 1명, c는 2명(정확히 절반)
    expr = pd.DataFrame({"a": [0, 2, 3, 1], "b": [0, 0, 0, 5], "c": [0, 0, 1, 1.5]})
    assert expressed_genes(expr, 0.5).tolist() == ["a", "c"]
