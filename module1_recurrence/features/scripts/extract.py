"""공통 wsi_pipeline(A 담당)을 호출해 TCGA 슬라이드를 feature로 바꾼다.

파일 이름과 형식은 A와 맞춘다. 모듈 2·3 외부 검증과 4-B에서도 같이 쓴다.
실행: python -m module1_recurrence.features.scripts.extract
결과: data/tcga/features/{patient_id}.h5
"""
