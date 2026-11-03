"""예측 결과 평가: Pearson·Spearman·MAE, 고위험 AUC, NPV ≥ 0.95 기준점, calibration.

주 평가는 recommended 코호트(N0 > 0.5cm, N1 폐경 후)로 따로 낸다.
실행: python -m module1_recurrence.scripts.evaluate
결과: outputs/module1/metrics.json, 그림
"""
