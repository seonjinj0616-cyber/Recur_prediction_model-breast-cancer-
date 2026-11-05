# 모듈 1 재발 검증용 보완 데이터 조사

작성일: 2026-10-30
관련 문서: `docs/research_datasets_models.md` (0. 데이터 선정 기준), `docs/module1_cohorts_plan.md` (모듈 1 계획 · 진행)

---

## 1. 왜 보완이 필요한가

| 항목 | TCGA-BRCA (현재) | 재발 검증에 필요한 수준 |
|---|---|---|
| 검사 권고 코호트 | 379명 | — |
| PFI 사건 | **25건** (high 그룹 4건) | 수십~수백 건 |
| 추적 중앙값 | 약 2년 | HR+ 유방암은 5년 이후 재발이 많아 **5–10년 이상** |
| 유사 RS의 예후 분리 | high 10.2% vs 나머지 7.0%, p = 0.21 (유의하지 않음) | — |
| 실제 Oncotype DX 점수 | 없음 | 있으면 가장 좋음 |

TCGA만으로는 "H&E → 유사 RS"까지는 검증할 수 있지만, **"그 점수가 실제 재발을 나눈다"**는 보일 수 없습니다.

---

## 2. 보완 전략: 두 단계로 나눠 검증

H&E와 재발 정보가 함께 있는 공개 데이터는 드뭅니다. 그래서 검증을 둘로 나눕니다.

```
[단계 A] 유사 RS 라벨 → 재발        : RNA + 장기 추적 공개 코호트 (METABRIC, SCAN-B)
[단계 B] H&E → 유사 RS               : TCGA (현재 계획)
[단계 C] H&E → 재발 (직접, 선택)     : H&E + 추적 코호트 (PLCO, ABCTB 등, 신청 필요)
```

A와 B가 모두 성립하면 "H&E 모델이 예후 정보를 가진 라벨을 재현한다"는 간접 근거가 됩니다. 직접 증거(C)는 신청형 데이터가 확보될 때 추가합니다.

---

## 3. 후보 데이터

### 3-1. 단계 A: RNA + 장기 추적 (H&E 없음, 바로 쓸 수 있음)
| 데이터 | 규모 | 내용 | 접근 | 쓰는 방법 |
|---|---|---|---|---|
| **METABRIC** | 약 1,900명 | 마이크로어레이 발현 + ER/HER2 + 재발 무생존(RFS) · 전체 생존, 추적 10년 이상 | cBioPortal(`brca_metabric`)에서 가입 없이 다운로드 | 같은 `pseudo_rs.py` 공식으로 RS 계산 → ER+/HER2- 환자에서 RS 그룹별 재발 곡선 · Cox 분석. **라벨 자체의 예후력 검증** |
| **SCAN-B** (GSE96058) | 약 3,000명 | RNA-seq + 임상, 전체 생존 | GEO 공개 | TCGA와 같은 RNA-seq라 척도가 비슷함. 생존 검증 보조 |

- H&E가 없어서 이미지 모델을 직접 검증하지는 못합니다.
- METABRIC은 마이크로어레이라 0–15 척도 맞추기를 METABRIC 안에서 따로 해야 합니다 (RS 순위 기준 그룹은 그대로 적용 가능).
- 규모·필드 수치는 기억에 기반한 값입니다. 다운로드 전에 cBioPortal / GEO 페이지에서 다시 확인하세요.

### 3-2. 단계 C: H&E + 추적 (신청 필요)
| 데이터 | 규모 | 결과 정보 | 접근 | 비고 |
|---|---|---|---|---|
| **PLCO Breast** (NCI CDAS) | 유방암 약 1,012명, H&E 2,909장 | 장기 사망·생존 추적 (검진 코호트라 추적이 김) | CDAS 프로젝트 신청 · 승인 | 승인된 프로젝트 중 H&E로 재발·후기 사망을 예측하는 연구가 여럿 있음. **가장 현실적인 직접 검증 후보** |
| **ABCTB** (Australian Breast Cancer Tissue Bank) | 수천 명 규모 | 재발·생존, 일부 Oncotype 점수 | 윤리·과학 심사 후 제공 | 2025 medRxiv 연구(Cohen 등)가 H&E → Oncotype 외부 검증에 사용 |
| **Nightingale – Providence** | 2,500명 이상 (생검 3,255건) | 암 등록·사망 자료 연결 | Nightingale 플랫폼 신청 | 생검 슬라이드라 모듈 2(BCNB 생검)와도 맞음. 데이터는 플랫폼 밖으로 반출 불가일 가능성 |

### 3-3. 공개지만 작거나 제약이 큰 데이터
| 데이터 | 규모 | 결과 정보 | 접근 | 한계 |
|---|---|---|---|---|
| **HistologyHSI-BC Recurrence** (TCIA, *Sci Data* 2025) | 47명, H&E WSI 47장 (+ 초분광 영상) | 원격 재발 22 / 무재발 25, 12년 추적, DFS · OS | TCIA 공개 | 너무 작고 ER+ 수가 문서에 없음. **작은 외부 sanity check** 용도 |
| CPTAC-BRCA (TCIA + PDC/GDC) | 약 134명, 642장 | 임상은 PDC/GDC, 재발 추적 여부 불확실 | 공개 | 추적 필드 확인 필요. 규모 작음 |
| Zenodo 6369837 | TCGA FFPE 88장 | 재발/전이 이진 라벨 5건 | 공개 (CC BY 4.0) | **TCGA에서 나온 것**이라 독립 검증 아님 |
| DRYAD ER+ (TCGA·CINJ·CWRU·HUP) | 584장 | 결과 정보 언급 없음 | 공개 (축소 이미지만) | 원본 크기 슬라이드 미공개 |
| Duke-Breast-Cancer-MRI (TCIA) | 922명 | **Oncotype 점수**, 재발 | 공개 | **MRI만 있고 H&E 없음** → 이미지 모델 검증 불가 |
| BCNB (확보 완료) | 1,058명 | ER/PR/HER2, 림프절 | 확보 | 재발·생존 추적 없음 |

---

## 4. 추천 순서

| 순위 | 할 일 | 이유 | 데이터 선정 기준(0절) |
|---|---|---|---|
| 1 | **METABRIC으로 유사 RS 라벨의 예후력 검증** | 가입 없이 바로 가능. 사건 수가 많아 통계적으로 의미 있는 결과를 낼 수 있음 | 공개 · 재배포 조건 확인 필요 |
| 2 | **HistologyHSI-BC로 작은 외부 확인** | 공개 H&E + 장기 재발 추적. 결과는 "참고용"으로만 보고 | TCIA 라이선스 확인 필요 |
| 3 | **PLCO CDAS 신청** (일정이 허락하면) | 독립 H&E + 장기 추적의 가장 현실적인 경로. 승인에 수 주~수 개월 | 승인 전까지 기준 1 미충족. 원본·라벨 재배포 금지일 가능성이 높아 BCNB처럼 코드·지표만 공개 |
| — | ABCTB, Nightingale | 실제 Oncotype·재발이 있어 이상적이지만 심사 부담이 큼 | 부트캠프 기간에는 비현실적 |

**보고서에 쓸 한계 문장 (예시)**
> 본 모듈은 TCGA-BRCA에서 RNA 기반 유사 RS를 H&E로 예측하도록 학습했다. 유사 RS의 예후력은 METABRIC에서 별도로 확인했으나, H&E 예측값과 재발의 직접 연관은 독립 장기 추적 코호트에서 검증하지 못했다.

---

## 5. 같은 방식의 선행 연구

| 연구 | 라벨 | 학습 | 외부 검증 | 코드 |
|---|---|---|---|---|
| **Howard et al., *npj Breast Cancer* 2023** | TCGA RNA로 계산한 research-based ODX. HR+/HER2- 안 상위 15%를 high로 둠 | TCGA 1,039명, Xception + 종양 영역 가중, Tennessee 임상 노모그램과 결합 | UCMC 427명 **실제 ODX** AUROC 0.83 (노모그램 0.76). 항암 미시행 322명 재발 C-index 0.743 (실제 ODX 0.776) | [fmhoward/DLRS](https://github.com/fmhoward/DLRS), 가중치 Zenodo 7490381. **라이선스 표기 없음** |
| Su et al., arXiv 2409.15491 (Deep-BCR-Auto) | Howard의 research RS 사용 | TCGA-BRCA | OSU 자체 데이터 AUROC 0.832 | 확인 못 함 |
| Moon et al., ASCO 2024 | TCGA 발현 기반 + 삼성서울병원 실제 점수 | TCGA 1,009 + SMC 483 | — | 확인 못 함 |
| Boehm et al. (Orpheus), *Nat Commun* 2025 | 실제 RS | 3개 기관 6,172건 | — | 확인 못 함 |
| Cohen et al., medRxiv 2025 | 실제 ODX | 5개 코호트 4,227명 | Sheba · ABCTB · UCMC AUC 0.82 | 확인 못 함 |

**우리 유사 RS와 Howard 라벨 비교** (DLRS 저장소의 `tcga_brca_complete.csv`, 겹치는 508명)
- 순위 상관 Spearman **0.86**
- Howard high 67명 중 59명이 우리 high, 8명이 intermediate. 우리 low 82명은 모두 Howard low

**재사용할 수 있는 것**
- Howard CSV에는 **조직학적 등급(Grade)과 종양 크기**가 있습니다. 우리 임상 파일에는 없던 값이라 임상 결합 모델에 쓸 수 있습니다 (크기는 구간값으로 보임).
- 저장소에 라이선스가 없어서 기본적으로 저작권이 유보된 상태입니다. 비교와 인용은 문제없지만, 코드나 CSV를 우리 공개 저장소에 그대로 넣으려면 저자에게 허락을 받아야 합니다.

## 6. 확인하지 못한 것
- METABRIC · SCAN-B의 정확한 환자 수와 필드는 원본 페이지에서 재확인하지 않았습니다.
- HistologyHSI-BC의 ER+ 환자 수와 데이터 라이선스, PLCO의 신청 소요 기간과 재배포 조건은 원문에서 확인하지 못했습니다.
- ABCTB 연구자 신청 페이지는 열지 못했습니다.
- ESMO 2023 초록(308P)은 "Oncotype 점수가 있는 TCGA 환자 100명"으로 외부 검증했다고 하지만, 그 점수의 공개 출처는 찾지 못했습니다.

## 출처
- [Publicly available datasets of breast histopathology H&E WSIs: scoping review (arXiv 2306.01546)](https://arxiv.org/html/2306.01546v2)
- [PLCO pathology images (CDAS)](https://cdas.cancer.gov/learn/plco/images/pathology) · [PLCO 승인 프로젝트 예: 유방암 예후](https://cdas.cancer.gov/approved-projects/3342)
- [Prediction of OncotypeDX recurrence score using H&E WSI (medRxiv 2025)](https://www.medrxiv.org/content/10.1101/2025.07.21.25331907)
- [ESMO 2023 Abstract 308P](https://oncologypro.esmo.org/meeting-resources/esmo-congress-2023/regression-based-deep-learning-predicts-breast-cancer-recurrence-risk-score-from-pathology-slides)
- [ABCTB 이미지 컬렉션 (researchdata.edu.au)](https://researchdata.edu.au/breast-cancer-tissue-image-collection/60182)
- [Nightingale – Identifying High-Risk Breast Cancer (Atlas)](https://atlaslongitudinaldatasets.ac.uk/datasets/identifying-high-risk-breast-cancer-using-digital-pathology-images-a-nightingale-open-science-dataset)
- [HistologyHSI-BC Recurrence (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12663337/)
- [CPTAC-BRCA (TCIA)](https://doi.org/10.7937/tcia.caem-ys80)
- [Zenodo 6369837](https://zenodo.org/records/6369837)
- [Duke Breast Cancer MRI dataset](https://sites.duke.edu/mazurowski/resources/breast-cancer-mri-dataset/)
