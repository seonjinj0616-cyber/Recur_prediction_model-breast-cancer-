"""환자 단위 5-fold + 별도 테스트 병원 묶음. 수집 병원(TSS)이 학습과 테스트에 겹치지 않게 한다.

실행: python -m module1_recurrence.splits
결과: data/tcga/splits.csv (patient_id, fold, is_test)
"""
