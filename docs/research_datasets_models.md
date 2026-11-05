# team_project_final — 데이터 선정 기준 · 모듈 1 배경

작성일: 2026-10-28 · 정리: 2026-11-05 (이 레포는 **모듈 1 재발 위험**만 다루므로 모듈 2 · 3 · 4-B, 팀 배정 · 일정, 사업화 아이디어는 뺐음. 이전 판은 git 기록에 있음)
모듈 1의 실제 계획 · 진행 · 결과는 `docs/module1_cohorts_plan.md`, 외부 코호트 조사는 `docs/recurrence_validation_data.md`

---

## 0. 데이터 선정 기준

이 프로젝트는 부트캠프 교육용으로 GitHub에 공개합니다. 그래서 아래 세 가지를 **모두** 만족하는 데이터만 씁니다.

1. 팀이 **실제로 접근 권한을 확보**했다
2. **라이선스를 원본 페이지에서 확인**했다
3. 라이선스가 **교육용 공개 저장소에 올릴 수 있는 범위**가 명확하다

### 0-1. 모듈 1이 쓰는 데이터
| 데이터 | 접근 방법 | 라이선스 | GitHub에 올릴 수 있는 것 |
|---|---|---|---|
| **TCGA-BRCA** (NCI GDC) | 가입 없이 다운로드 (open access) | NIH GDC open-access 데이터 (인용 필요) | 슬라이드 ID 목록, 라벨 CSV, 추출한 feature |
| **TCGA-CDR** (Liu et al., *Cell* 2018 부록) | 논문 부록 다운로드 | 논문 부록 공개 데이터 (인용 필요) | 사용한 환자 ID와 endpoint 열 |
| METABRIC · SCAN-B · SEQUOIA 서명 | `docs/module1_cohorts_plan.md` 1-5 참고 | 데이터마다 다름 (같은 표 참고) | 같은 표 참고 |

- **원본 WSI는 GitHub에 올리지 않습니다.** 파일 하나가 수백 MB~수 GB라 GitHub 100MB 제한에 걸립니다.
- feature 추출 encoder는 **torchvision ImageNet ResNet50**을 씁니다. 병리 전용 모델(UNI, CONCH 등)은 접근 승인과 비상업 약관이 붙어서 기본값에서 뺐습니다.
- 세포 분할(CellViT-SAM-H)은 PanNuke 학습 가중치라 **비상업 연구용**입니다 (`module1_cohorts_plan.md` 단계 1).

### 0-2. GitHub 공개 원칙
| 올리는 것 | 올리지 않는 것 |
|---|---|
| 다운로드 스크립트, 슬라이드 ID manifest (TCGA) | 원본 WSI (용량 초과) |
| 라벨 CSV (TCGA 임상에서 뽑은 열) | 라이선스상 재배포가 안 되는 데이터 (METABRIC · Rueda 부록 등, 계획 문서 1-5 표) |
| feature `.h5` (TCGA. 용량이 크면 Hugging Face 등 별도 저장소) | 병리 foundation model 가중치 |
| 학습 코드, 학습된 모델 가중치 | API 키, 개인 정보 |

README에 데이터별 **출처, 라이선스, 인용 문구**를 적습니다.

---

## 1. 프로젝트 안에서 모듈 1의 위치

| # | 모듈 | 입력 → 출력 |
|---|---|---|
| **1** | **재발 위험 (이 레포)** | H&E → 원격 재발 위험, Oncotype 검사 권고 근거 |
| 2 | 수술 전 림프절 전이 예측 | 생검 H&E (+임상) → 림프절 전이 확률 |
| 3 | 바이오마커 예측 | H&E → ER · PR · HER2 · Ki-67 |
| 4-B | 멀티모달 유사 환자 검색 | H&E 슬라이드 벡터 + 임상 → 비슷한 과거 환자 |

**전체 흐름:** 생검 H&E → 모듈 1 · 2 · 3 → 모듈 4-B → 환자 종합 리포트. 슬라이드 전처리(`wsi_pipeline`)와 feature는 조원마다 각자 만듭니다.

---

## 2. 용어

| 용어 | 쉬운 설명 |
|---|---|
| **H&E 사진** | 조직을 얇게 썰어 분홍색(E, 에오신)과 보라색(H, 헤마톡실린)으로 물들인 현미경 사진. 암 진단의 기본 |
| **WSI (Whole Slide Image)** | 유리 슬라이드 한 장 전체를 고해상도로 스캔한 사진. 수만 × 수만 픽셀이라 작은 조각(**패치 · 타일**)으로 잘라서 씀 |
| **IHC 검사** | 특정 단백질만 갈색(DAB)으로 물들이는 특수 염색. TCGA 진단 슬라이드 중 1장이 H&E 대신 IHC였음 (계획 문서 슬라이드 QC) |
| **Oncotype DX** | 암 조직의 유전자 21개로 재발 위험 점수(0–100)를 내는 검사. 이 점수로 항암치료 여부를 정함. 매우 비쌈 |
| **유사 RS** | 공개 데이터에 실제 Oncotype 점수가 없어서 RNA-seq에 공개된 21-gene 공식을 적용해 직접 계산한 점수 |
| **MIL (Multiple Instance Learning)** | 타일마다 정답이 없어도 "환자 하나의 정답"만으로 학습하는 방법. 어떤 타일이 중요했는지 히트맵으로 보여줄 수 있음 |

---

## 3. 기존 모델 (H&E → Oncotype · 재발)

- **Orpheus** (Boehm et al., *Nat Commun* 2025): H&E + 병리 리포트 멀티모달 transformer. RS ≤ 25 환자에서 원격 재발 예측 td-AUC 0.75로 RS 자체(0.49)보다 높음.
- **Deep-BCR-Auto**: attention MIL, TCGA AUROC 0.827.
- **AnchorMIL** (2025): 회귀 + 분류 동시 학습, TCGA AUC 0.89.
- Prov-GigaPath fine-tuning (Cohen et al., medRxiv 2025).
- **SEQUOIA** (2024): H&E → 유전자 발현 예측 → 재발 서명. 모듈 1의 가장 직접적인 비교 대상 (계획 문서 2장).
- 2026 체계적 문헌고찰: 연구 11편이 있지만 무작위 임상시험 데이터로 검증된 것은 없음.

**한계:** 유사 RS는 진짜 RS가 아님. TCGA는 추적이 짧고(약 2년) 재발 사건이 적어 C-index가 흔들릴 수 있음 → 원격 재발은 METABRIC · SCAN-B로 봄 (계획 문서).

---

## 4. 참고 링크
- **TCGA:** [NCI GDC Data Portal](https://portal.gdc.cancer.gov/) · TCGA-CDR: Liu et al., *Cell* 2018 (TCGA Pan-Cancer Clinical Data Resource)
- **Oncotype 예측 모델:** [Orpheus (Nat Commun 2025)](https://pubmed.ncbi.nlm.nih.gov/40025017/) · [AnchorMIL](https://research.itu.edu.tr/en/publications/hampe-to-oncotype-dx-predicting-recurrence-risk-in-hrher2-breast-/) · [Deep-BCR-Auto](https://aitopics.org/tag/tcga-brca%20dataset) · [Cohen et al. medRxiv 2025](https://www.medrxiv.org/content/10.1101/2025.07.21.25331907) · [2026 리뷰](https://makadu.live/packages/uploads/2026/03/2026-Deep-learning-AI-for-prognosis-prediction.pdf)
- **H&E 기반 분자 예측 리뷰:** [arXiv 2412.10392](https://arxiv.org/pdf/2412.10392) · [공개 유방 H&E WSI 데이터셋 리뷰 (PMC10884505)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10884505/)
