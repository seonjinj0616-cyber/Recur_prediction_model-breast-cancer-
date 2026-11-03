"""pseudo_rs 공식 확인: 손으로 계산한 예시 값, 0–100 자르기, 유전자 이름 변환(STK15→AURKA 등)."""
import numpy as np
import pandas as pd
import pytest

from module1_recurrence.labels.scripts.pseudo_rs import (CANCER_GENES, REFERENCE, assign_risk,
                                                 normalize, score)


def flat(value, **override):
    s = pd.DataFrame({g: [float(value)] for g in CANCER_GENES})
    for g, v in override.items():
        s[g] = float(v)
    return s


def test_hand_calculated():
    # 모두 9: 그룹 점수 모두 9, RSU = 9 × (0.47 − 0.34 + 1.04 + 0.10 + 0.05 − 0.08 − 0.07) = 10.53
    r = score(flat(9)).iloc[0]
    assert r.rs_unscaled == pytest.approx(10.53)
    assert r.pseudo_rs == pytest.approx(20 * (10.53 - 6.7))


def test_floors():
    # 모두 0: HER2 그룹 8, 증식 그룹 6.5로 올라감 → RSU = 0.47·8 + 1.04·6.5 = 10.52
    r = score(flat(0)).iloc[0]
    assert (r.her2_gs, r.proliferation_gs) == (8, 6.5)
    assert r.rs_unscaled == pytest.approx(10.52)


def test_clip_0_100():
    high = score(flat(15, ESR1=0, PGR=0, BCL2=0, SCUBE2=0)).pseudo_rs.iloc[0]
    low = score(flat(0, ESR1=15, PGR=15, BCL2=15, SCUBE2=15)).pseudo_rs.iloc[0]
    assert (high, low) == (100, 0)


def test_aliases():
    old = flat(9).rename(columns={"AURKA": "STK15", "CTSV": "CTSL2", "BIRC5": "Survivin"})
    pd.testing.assert_frame_equal(score(old), score(flat(9)))


def test_normalize_range_and_reference():
    rng = np.random.default_rng(0)
    log_expr = pd.DataFrame(rng.normal(5, 1, (20, 21)), columns=CANCER_GENES + REFERENCE)
    s = normalize(log_expr)
    assert s.min().eq(0).all() and s.max().eq(15).all()
    # 기준 유전자가 한 샘플에서 모두 같이 오르면(전체 발현량 차이) 정규화 값은 그대로
    shifted = log_expr.copy()
    shifted.iloc[0] += 3
    pd.testing.assert_frame_equal(normalize(shifted), s)


def test_assign_risk_proportions():
    groups = assign_risk(pd.Series(np.arange(100.0)))
    assert groups.value_counts().to_dict() == {"intermediate": 67, "low": 16, "high": 17}
    assert groups.iloc[-1] == "high" and groups.iloc[0] == "low"
