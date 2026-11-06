# 모듈 1 (H&E 기반 원격 재발 위험) — 코호트 · 비교 모델 · 개발 단계

작성일: 2026-10-30 · 수정: 2026-11-05 (진행 정리 문서 `module1_progress_summary.md`를 이 문서로 합침: 6장 처리 규칙, 7장 실행 순서, 8장 용어)
관련 문서: `docs/recurrence_validation_data.md` (외부 코호트 조사), `docs/research_datasets_models.md` (데이터 선정 기준 · 배경)

모듈 목표: **H&E에서 읽을 수 있으면서 원격 재발과 관련된 유전자**로 재발 위험을 계산하고, Oncotype DX 검사 대상 중 **임상적으로 놓치기 쉬운 고위험 환자에게 검사를 권고하는 근거**를 준다. (Oncotype 대체가 아님)

---

## 1. 데이터별 코호트

공통 규칙: HR+ = ER 또는 PR 양성, HER2 음성. 검사 권고 = **N0 & 종양 > 0.5 cm** 또는 **N1(1–3개) & 폐경 후** (pN1mi는 N0).

### 1-1. TCGA-BRCA — H&E ↔ RNA 연관 (단계 1·2), 이후 H&E 평가
| 코호트 | 환자 수 | 용도 |
|---|---|---|
| 여성 + 진단 슬라이드 + RNA-seq (아형 무관) | **1,047** (HR+/HER2− 531, 그 외 516) | 단계 2: H&E → 유전자 연관 (SEQUOIA와 같은 조건: 유방암 전체) |
| HR+/HER2− (M0, 선행항암 없음, 슬라이드 + RNA) | **515** | 단계 2 민감도 분석, 이후 H&E 재발 평가 |
| ↳ 검사 권고 (N0 278 + N1 폐경 후 101) | **379** | H&E 단계 주 평가 |

| HR+/HER2− 그룹 | n | PFI 사건 | 5년 내 | 추적 중앙값 |
|---|---|---|---|---|
| N0 | 278 | 14 | 11 | 2.2년 |
| N1 폐경 후 | 101 | 11 | 11 | 1.6년 |
| N1 폐경 전 | 52 | 5 | 4 | 2.0년 |
| N2+ | 80 | 8 | 7 | 1.9년 |
| 림프절 미상 | 4 | 1 | 1 | 2.1년 |
| **합계** | **515** | **39** | **34** | 2.0년 |

- 재발 정의: TCGA-CDR PFI (국소 + 원격 + 진행). **원격 재발만 따로 없음.**
- 현재 RNA-seq은 HR+/HER2− 515명분만 받음 → 단계 2를 1,047명으로 하려면 나머지 532명 RNA와 전체 슬라이드가 필요.

### 1-2. METABRIC — RNA ↔ 원격 재발 (단계 3), 서명 학습
Rueda et al. 2019 부록(Table S5)의 **원격 재발(DR, TDR)** 사용. HR+/HER2− 1,413명 모두 mRNA(Illumina HT-12 마이크로어레이) 있음.

| 그룹 | 전체 n | 원격 재발 | 10년 내 | 추적 중앙값 | 재발 기록 완전 n | 원격 재발 |
|---|---|---|---|---|---|---|
| N0 (> 0.5 cm) | 742 | 154 | 127 | 10.5년 | 664 | 102 |
| N1 폐경 후 | 368 | 112 | 98 | 8.6년 | 325 | 77 |
| N1 폐경 전 | 56 | 16 | 11 | 11.1년 | 51 | 11 |
| N2+ | 194 | 94 | 83 | 6.1년 | 169 | 72 |
| T ≤ 0.5 cm | 13 | 4 | 3 | — | 7 | 1 |
| 림프절 미상 | 40 | 12 | 9 | — | 35 | 9 |
| **HR+/HER2− 합계** | **1,413** | **392** | **331** | 9.5년 | **1,251** | **272** |

| 하위 코호트 | n | 원격 재발 | 10년 내 |
|---|---|---|---|
| 검사 권고 | 1,110 | 266 | 225 |
| 검사 권고 + 재발 기록 완전 (**주 분석**) | **989** | **179** | **154** |
| ↳ 내분비요법만 (Oncotype 개발 조건과 비슷) | 632 | 108 | 94 |

- 치료 (검사 권고 + 기록 완전): 내분비만 632, 무치료 298, 항암 + 내분비 51, 항암만 8.
- 그룹은 `cohort.py`의 `assign_group` 규칙(TCGA와 같음)으로 나눔: 종양 ≤ 0.5 cm이면 림프절 상태와 상관없이 T ≤ 0.5 cm 그룹. 종양 크기 결측은 > 0.5 cm로 처리.
- 다른 세션 집계(검사 권고 1,102명, 원격 재발 264건)와 약간 다른 이유: 종양 크기 결측 처리 차이.
- 표 파일: `data/metabric/cohort.csv` (`python -m module1_recurrence.data.scripts.external_cohort`)

### 1-3. SCAN-B (GSE96058) — 서명 외부 검증 (단계 5)
RNA-seq, 원격 재발 무사건 구간(**DRFi**) 있음. 추적 코호트(Follow.up.cohort) 환자당 1샘플.

| 그룹 | n | 원격 재발 | 5년 내 | 추적 중앙값 |
|---|---|---|---|---|
| N0 (> 0.5 cm) | 2,340 | 85 | 60 | 5.3년 |
| N1 "폐경 후" | 755 | 43 | 28 | 5.3년 |
| N1 "폐경 전" | 231 | 16 | 10 | 5.5년 |
| N2+ | 324 | 75 | 63 | 5.2년 |
| T ≤ 0.5 cm / 미상 | 114 | 3 | 2 | — |
| **HR+/HER2− (DRFi 있음)** | **3,764** | **222** | **163** | 5.3년 |
| **검사 권고** | **3,095** | **128** | — | — |
| ↳ 내분비요법만 | 2,076 | 81 | — | — |

- 표 파일: `data/scanb/cohort.csv` (HR+/HER2− 추적 코호트 전체, DRFi 없는 환자 포함)

- **폐경 정보가 없어서 나이 구간 51세 이상을 폐경 후로 봄** (근사).
- PAM50 ROR(NCN.ROR, SSP.ROR)이 이미 계산되어 있어 비교 기준으로 바로 쓸 수 있음.

### 1-4. CPTAC-BRCA — H&E → 유전자 외부 검증 (단계 2 확인용)
| 항목 | 내용 |
|---|---|
| RNA-seq | 133명, GDC 공개, TCGA와 같은 STAR-Counts 형식 |
| 슬라이드 | TCIA 642장 (약 113 GB) |
| ER/HER2 | GDC에 없음 → Krug 2020 Table S1에서 가져와야 함 (미확인) |
| 재발 정보 | 없음 → 재발 검증 불가 |

### 1-5. 데이터 위치
다운로드: `python -m module1_recurrence.data.scripts.external_fetch` (발현 행렬까지: `--expression`, 주소만 확인: `--check`)

| 경로 | 내용 | 출처 · 라이선스 | GitHub |
|---|---|---|---|
| `data/metabric/raw/clinical_patient.json`, `clinical_sample.json`, `mrna_sample_ids.json` | 임상, mRNA 있는 샘플 1,980개 | cBioPortal `brca_metabric` | 올리지 않음 |
| `data/metabric/raw/rueda2019_TableS5.txt`, `TableS6.txt` | 원격 재발 DR/TDR, 국소 재발, 재발 기록 완전 여부 | Rueda et al. *Nature* 2019 부록 (라이선스 표기 없음) | 올리지 않음 |
| `data/metabric/hr_her2neg.csv` | 위 임상에서 만든 HR+/HER2− 1,413명 표 | 파생 | 올리지 않음 |
| `data/scanb/raw/scanb_clinical_2023-01-13.xlsx` | 9,206 샘플 임상 + DRFi + PAM50 ROR | Staaf et al. 2022, Mendeley Data 10.17632/yzxtxn4nmd.3, **CC BY 4.0** | 출처 표기 시 가능 |
| `data/scanb/raw/Gene.ID.ann.txt` | 유전자 주석 (GENCODE 27) | 같은 Mendeley, CC BY 4.0 | 출처 표기 시 가능 |
| `data/scanb/raw/GSE96058-*_series_matrix.txt.gz` | GEO 샘플 메타데이터 | GEO GSE96058 | 가능 |
| `data/metabric/raw/data_mrna_illumina_microarray.txt` | 20,603행 × 1,980샘플 (Illumina HT-12) | cBioPortal datahub | 올리지 않음 (용량) |
| `data/scanb/raw/SCANB.9206.genematrix_noNeg.txt` | 19,675유전자(ENSG) × 9,206샘플, 라이브러리 보정 FPKM | 같은 Mendeley, CC BY 4.0 | 올리지 않음 (용량) |

- HR+/HER2− 코호트 전원이 발현 행렬에 있음: METABRIC 1,413/1,413, SCAN-B(DRFi 있음) 3,764/3,764.
- Oncotype 21개 유전자는 두 행렬 모두에 있음 (CTSL2는 현재 이름 **CTSV**로 찾음).
- METABRIC은 같은 유전자 기호가 여러 행(프로브)인 경우가 219개 → 유전자별로 평균내서 사용.
- SCAN-B 행은 Ensembl ID → `Gene.ID.ann.txt`로 유전자 기호 변환.
- 전처리 결과: `data/metabric/expr_hr_her2neg.parquet`, `data/scanb/expr_hr_her2neg.parquet`, 공통 유전자 목록 `data/expr_common_genes.txt` (`python -m module1_recurrence.data.scripts.expression`)

---

## 2. 비교할 기존 모델

### 2-1. RNA 단계 (METABRIC 교차검증 · SCAN-B에서 같은 환자로 직접 비교)
| 비교 대상 | 기능 | 계산 방법 | 답하는 질문 |
|---|---|---|---|
| **유사 RS** (Oncotype 공식) | 유전자 검사 | `pseudo_rs.py` 공식 (마이크로어레이는 코호트 안에서 척도 맞춤) | Oncotype 공식보다 원격 재발을 잘 나누나 |
| **PAM50 ROR** (Prosigna) | 유전자 검사 | SCAN-B 제공값, METABRIC은 genefu로 계산 | 다른 상용 서명과 비교 |
| **SEQUOIA 272개 서명** | H&E-예측 가능 유전자로 만든 재발 서명 (TCGA 학습) | 논문 Supplementary Data 5의 계수를 그대로 적용 | **가장 직접적인 경쟁 상대**: 재발 유전자를 TCGA가 아닌 METABRIC에서 고른 효과 |
| NPI, CTS5 | 임상 도구 | METABRIC NPI 열, CTS5는 크기·림프절·등급·나이로 계산 | 임상 정보만으로 충분한가 |
| 증식 점수, 무작위 유전자 세트 | 확인용 기준선 | Oncotype 증식 5개 평균 / 같은 개수 무작위 × 1,000 | 증식만 배운 것 아닌가 |

### 2-2. H&E 단계 (TCGA, 같은 fold · 같은 특징으로 재학습해서 비교)
| 비교 대상 | 기능 | 출처 |
|---|---|---|
| H&E → 재발 직접 학습 | H&E → 재발 위험 | SEQUOIA 대조군 방식 |
| H&E → 유사 RS | H&E → 유전자 점수 | Howard 2023 방식 |
| H&E → 유전자 → SEQUOIA 서명 | H&E → 유전자 점수 | SEQUOIA 2024 방식 |
| AI 유사분열 수, TIL 밀도 | H&E → 형태 지표 | 세포 분할 결과에서 계산 |

### 2-3. 문헌 수치만 참고 (직접 비교 불가: 비공개 또는 다른 데이터)
| 모델 | 입력 | 종점 | 보고 성능 |
|---|---|---|---|
| Orpheus (Boehm, *Nat Commun* 2025) | H&E + 병리 리포트 | 원격 재발 | td-AUC 0.75 (RS ≤ 25 환자) |
| RlapsRisk BC (Owkin) | H&E + 임상 | 5년 무전이 생존 | C-index 0.80 (임상만 0.76) |
| Ataraxis Breast | H&E (+ 임상) | 무병 기간 | C-index 0.71 (회사 발표) |
| DeepGrade / Stratipath | H&E | 무재발 생존 · 원격 재발 | C-index 0.71 |
| TCGA-BRCA MIL 벤치마크 | H&E | 전체 생존 위주 | C-index 0.54–0.64 |

수치는 원문 확인 수준이 제각각이므로 발표 전에 원문으로 다시 확인합니다.

---

## 3. 개발 단계 (1–5)와 기존 논문과의 차별점

| 단계 | 내용 | 데이터 · 코호트 | 기존 논문과 다른 점 |
|---|---|---|---|
| **1. 종양 밀집 영역만 사용** | 세포 분할 모델로 타일별 종양 세포 비율을 계산하고, 종양 비율이 높은 타일만 남김 | TCGA 1,047명 슬라이드 | SEQUOIA는 배경만 빼고 **무작위 타일** 사용, 종양 선별 없음 |
| **2. H&E ↔ RNA 연관 순위** | 남긴 타일로 유전자 발현 예측 → **SEQUOIA와 같은 기준**으로 "잘 예측된 유전자" 판정 → Pearson r 순으로 순위 | TCGA 1,047명 (HR+/HER2− 531명은 민감도), CPTAC로 외부 확인 | 기준은 SEQUOIA와 같게 두고, **순위 기준(Pearson r)을 명시** (SEQUOIA는 미기재) |
| **3. RNA ↔ 원격 재발** | 유전자별 Cox (림프절·크기·등급·나이 보정) → FDR < 0.05 전부를 R 목록으로 | METABRIC HR+/HER2− 1,251명 (기록 완전, 원격 재발 272건) | SEQUOIA는 재발 유전자를 **TCGA(추적 약 2년, 재발 전체)**에서 고름 → 우리는 **장기 추적 + 원격 재발** |
| **4. 교집합** | P(단계 2) ∩ R(단계 3), 두 플랫폼에 모두 있는 유전자만 | — | 두 목록을 **독립 코호트에서** 만든 뒤 교집합 (SEQUOIA는 같은 TCGA에서 순차 필터) |
| **5. 재발 수치 계산** | 교집합 유전자로 Cox 모델 학습 → 가중치 고정 → 재발 위험 점수 + 10년 원격 재발 확률 → SCAN-B 외부 검증 | 학습: METABRIC / 검증: SCAN-B 검사 권고 3,095명 (원격 재발 128건) | HR+/HER2−·검사 권고 코호트 평가, 유사 RS와 같은 환자로 비교 |

### 단계별 세부 기준

**단계 1. 종양 영역 선별**
- 도구: **CellViT++의 CellViT-SAM-H** (2026-11-04 선정, TCGA 슬라이드 시험 후 확정). PanNuke 5종 (종양 · 상피 · 염증 · 결합조직 · 괴사) 중 **종양(Neoplastic)** 핵 비율을 씀.
- 기준 예시: 타일 안 종양 핵 비율 ≥ 50% 타일만 사용 (기준값은 분포를 보고 결정, 결정 후 고정).
- 주의: 면역·기질 유전자는 덜 예측될 수 있음 → **"모든 조직 타일" 버전도 함께 돌려** 교집합이 어떻게 바뀌는지 비교.
- 분할은 픽셀이 필요해서 `.h5` 특징만으로는 다시 할 수 없음 → 슬라이드마다 특징 추출과 같이 돌려 **타일 좌표 · 특징 · 종양 핵 비율을 한 `.h5`에** 저장한 뒤 svs 삭제. 기준값을 바꿔도 다시 받을 필요 없음.

**세포 분할 모델 선정 (단계 1)**

선정 기준: ① 어느 핵을 종양으로 판정했는지 눈으로 확인 가능 (WSI 처리 + QuPath 등에서 보기) ② 공개된 모델 중 성능 ③ 최근 발표 ④ 읽을 수 있는 논문.

| 모델 | ① 확인 | PanNuke mPQ / bPQ / 종양 PQ | 발표 | 코드 · 가중치 |
|---|---|---|---|---|
| **CellViT++ (SAM-H)** | `.svs` 직접, GeoJSON, 웹 뷰어 | 0.498 / 0.679 / 0.581 | 2025.01 (arXiv 2501.05269) | 공개, 비상업 |
| CellViT-SAM-H | GeoJSON | 0.498 / 0.679 / 0.581 | 2024 (MedIA 94, 103143) | 공개, 비상업 |
| NuLite-H | WSI, GeoJSON | 0.496 / 0.677 / 0.577 | 2024.08 (BSPC 2026) | 공개, 비상업 |
| LKCell-L | 없음 (PanNuke 이미지용만) | 0.508 / 0.685 / 0.586 | 2024.07 (Neurocomputing) | 공개 (MIT) |
| CFR-SAM-H | 코드 없음 | 0.511 / 0.696 / 0.599 | 2026.03 (MedIA) | 확인 못 함 |
| NucShare-T | 없음 | 0.524 / 0.681 / — | 2026.09 프리프린트 | 가중치 요청 필요 |
| HoVer-Net | 타일만 | 0.463 / 0.660 / 0.551 | 2019 (MedIA) | PyTorch 1.6이라 RTX 5090 불가 |

- PQ = 탐지 F1 × 맞힌 핵의 평균 IoU (0–1). mPQ는 5종 평균, bPQ는 종류 무시. 수치 출처: CellViT++ 부록 Table S1·S2, LKCell Table 2·3, NuLite Table 2·4, CFR-SAM · NucShare Table 1.
- **선정 이유:** 비상업 라이선스 모델 중 **코드와 가중치가 오픈소스로 공개**되어 있고 ①을 만족함. CellViT++ 논문에서 비교한 인코더 5개 (SAM-H, HIPT256, UNI, Virchow, Virchow2) 중 **SAM-H가 PanNuke mPQ가 가장 높음** (0.498, 나머지 0.485–0.493).
- 더 높은 CFR-SAM · NucShare는 코드나 가중치를 바로 쓸 수 없고, LKCell은 WSI 처리 · 시각화가 없음. NuLite-H는 성능이 비슷하고 더 가벼워서 처리 시간이 문제가 되면 대안.
- 라이선스: PanNuke가 CC BY-NC-SA 4.0이라 가중치는 비상업 연구용. CellViT++ 코드는 Apache 2.0 + Commons Clause (판매 금지, 인용 필수). README상 일부 체크포인트(UNI · Virchow 계열로 추정)는 공유되지 않지만 SAM-H는 Google Drive와 Docker 이미지로 제공.
- 요구 사항: GPU 24 GB 이상 (RTX 5090 32 GB 충족), Python 3.10.
- 실행: PyPI `cellvit` 패키지 (`cellvit-inference --model SAM --nuclei_taxonomy pannuke --geojson process_wsi --wsi_path …`). 공식 안내의 PyTorch 2.2(cu121)는 RTX 5090을 지원하지 않아 **torch 2.11 + cu128**로 설치, 별도 conda 환경 `cellvit`. 가중치(2.8 GB)는 첫 실행 때 `~/.cache/cellvit`에 받음.

**CellViT 시험 결과 (2026-11-04, GPU 서버, QuPath로 눈 확인)**

| 슬라이드 | 염색 · 특징 | 패치 (1024) | 시간 | 핵 수 (Neoplastic 비율) | 판정 |
|---|---|---|---|---|---|
| TCGA-5L-AAT0 | **IHC** (상피 표지자 DAB, 에오신 없음) | 325 | 6.5분 | 105,239 (25%) | 기질까지 종양으로 판정 → 사용 불가. **염색 문제** |
| TCGA-A1-A0SB | H&E, 초록 · 파란 펜 | 653 | 15분 | 381,686 (50%) | 종양 덩어리 · 기질 · 정상 상피 구분 맞음. **진한 펜 아래 핵 놓침, 옅은 파란 펜 아래 기질 핵을 종양으로 오분류** |
| TCGA-5T-A9QA | H&E, 깨끗함 | 1,717 | 55분 (후처리 12분), 메모리 최대 45 GB | 1,814,963 (58%) | 종양 · 기질 · 림프구 구분 맞음, 고배율에서 핵 대부분 검출 |

- 결론: H&E에서는 쓸 만함. **IHC 슬라이드 제외**와 **펜 영역 타일 제외**가 필요 → `data/scripts/slide_qc.py` (GDC 타일 API 썸네일로 염색 · 펜 판정).
- 처리 시간은 패치 수보다 빠르게 늘어남 (후처리). 515장은 GPU 2장으로도 며칠 단위, 큰 슬라이드는 메모리 45 GB라 동시 실행 개수 확인 필요.
- 초록 펜은 대부분 종양 위치 표시 → 선택 기준이 아니라 **판정 확인용 보조 근거**로만 사용.

**종양 범위 지정: CellViT vs BCSS 비교 (2026-11-05) → CellViT 사용, BCSS 사용 안 함**

비교 코드(`features/scripts/bcss_compare.py`, 결정 후 삭제 · git 기록에 있음): TIAToolbox `fcn_resnet50_unet-bcss` (BCSS로 학습한 조직 영역 분할, 가중치 CC0)를 같은 128 µm 타일 격자에서 돌려
(256 µm 창 → 0.25 mpp 1024픽셀, 가운데 512픽셀 출력 = 타일 하나) 타일별 Tumour 면적 비율 ≥ 0.5를 BCSS 종양 타일로 보고, CellViT 종양 타일(종양 핵 비율 ≥ 0.5, 핵 ≥ 20)과 비교.

| 슬라이드 | 특징 | CellViT 종양 | BCSS 종양 | 일치율 | kappa | Dice |
|---|---|---|---|---|---|---|
| TCGA-5T-A9QA | 깨끗함, 종양 덩어리 | 12,680 | 10,454 | 0.85 | 0.69 | 0.87 |
| TCGA-A1-A0SB | 펜, 종양 덩어리 | 1,374 | 968 | 0.84 | 0.64 | 0.76 |
| TCGA-A8-A0AD | 바랜 H&E | 501 | 752 | 0.80 | 0.60 | 0.77 |
| TCGA-A8-A096 | 가장 바랜 H&E | 437 | 337 | 0.77 | 0.49 | 0.66 |
| TCGA-A1-A0SD | 기질 속에 퍼진 종양 | 653 | 171 | 0.85 | 0.30 | 0.35 |
| TCGA-A1-A0SE | 기질 속에 퍼진 큰 종양 (초록 펜 표시) | 2,191 | 239 | 0.81 | 0.16 | 0.19 |

- 종양이 덩어리진 슬라이드에서는 두 방법이 잘 맞음 (Dice 0.76–0.87). CellViT가 종양 안 옅은 부위 · 경계를 조금 더 넓게 잡음.
- 기질 속에 퍼진 종양에서는 BCSS가 종양을 거의 놓침: A1-A0SE는 병리의가 펜으로 표시한 큰 종양을 BCSS가 239타일만 종양으로 봄.
- 원인은 기준 차이: BCSS는 **면적** (CellViT 종양 타일 안을 BCSS는 기질 63% · 종양 30%로 봄), CellViT는 **핵 비율** (같은 타일에서 종양 핵 68%). 기질은 핵이 적어 핵 비율로는 종양이 드러남.
- **결정:** CellViT 기준(당시 종양 핵 비율 ≥ 0.5, 핵 ≥ 20. 이후 6-1에서 RNA 조각 맞춤으로 바꿈)을 그대로 쓰고 BCSS는 쓰지 않음. 비교 결과 파일 · BCSS 환경 · 비교 코드(`bcss_compare.py`)는 지움. 코드는 git 기록에 있음.

**특징 인코더 (단계 2 입력) → H-optimus-0 하나로 확정 (2026-11-06)**

라이선스 기준(승인 없이 받기, 특징 공개 가능)을 지키는 후보 중 병리 전용 2개와 ImageNet ResNet50을 시험 실행으로 비교해 H-optimus-0 하나로 정함. ResNet50 · Hibou-L 코드와 시험 코드(`encoder_pilot.py`)는 지움 (git 기록에 있음).

| 모델 | 구조 · 차원 | 학습 데이터 | 라이선스 · 받기 | 벤치마크 (eva 9과제 평균) | 판단 |
|---|---|---|---|---|---|
| ImageNet ResNet50 | layer3 + 평균 풀링, 1024 | ImageNet | BSD, torchvision | — | 시험에서 탈락 (유전자 과제가 무작위와 같음) |
| **H-optimus-0** (Bioptimus) | ViT-g (약 11억), 1536 | 자체 50만 장 이상 | Apache 2.0, HF 로그인 + 연락처 동의 | 0.793 | **채택** (`extract.py`) |
| **Hibou-L** (HistAI) | ViT-L | 자체 100만 장 이상 | Apache 2.0 (`cellvit-hibou-l` 등 파생은 비상업) | 0.777 | 시험에서 탈락 |
| Phikon / Phikon-v2 (Owkin) | ViT-B / ViT-L | **TCGA** | Owkin 비상업 | 0.753 | 제외: 비상업 + TCGA 평가 시 누수 |
| CTransPath | Swin-T, 768 | **TCGA** + PAIP | GPL-3.0 + 비상업 학술 문구 | 최신 모델보다 낮음 | 제외: 같은 누수 · 옛 `timm` 필요 |

- 벤치마크 출처: MINT (arXiv 2603.07895) eva 표, Neidlinger et al. *Nat Biomed Eng* 2025 (19개 모델, H-optimus-0 상위권). 라이선스는 2026-11-05 Hugging Face 모델 페이지 기준. H-optimus-0 · Hibou-L 모두 HF 로그인 + 이용 동의 필요 (자동 승인).
- **시험 방법 (2026-11-05):** 처리된 50장에서 RNA 조각 맞춤 타일 · 그 밖의 조직 타일을 슬라이드당 300개씩 (27,433타일) 세 인코더에 넣고 ① GPU 속도 ② 타일 미니 과제: CellViT 핵 구성을 ridge로 맞힘, 슬라이드 단위 5-fold R² ③ 실제 과제 축소판: 환자별 맞춤 타일 평균 → R 목록 493개 · Oncotype 16개 log2(TPM+1), 환자 5-fold ridge, 유전자별 Pearson r. 교차검증 r은 신호가 없어도 음수로 나오므로 무작위 특징과 비교.

| | ResNet50 | **H-optimus-0** | Hibou-L | 무작위 (기준선) |
|---|---|---|---|---|
| 차원 | 1024 | 1536 | 1024 | 1024 |
| GPU 속도 (타일/초) | 6,227 | **224** | 708 | — |
| 510장 예상 (조직 타일 평균 7,693개 × 510) | 약 1.5시간 (읽기 한계) | **약 5시간** | 약 1.7시간 | — |
| 종양 핵 비율 R² | 0.814 | **0.868** | 0.829 | −0.003 |
| 염증 세포 비율 R² | 0.695 | **0.819** | 0.611 | −0.003 |
| 기질 세포 비율 R² | 0.783 | 0.838 | **0.850** | −0.004 |
| log 핵 수 R² | 0.897 | 0.930 | **0.942** | −0.002 |
| 유전자 발현 r 중앙값, R 목록 493개 | −0.314 | **0.120** | −0.173 | −0.320 |
| r > 0.3 유전자 수, R 목록 | 0 | **82** | 0 | 0 |
| 유전자 발현 r 중앙값, Oncotype 16개 | −0.326 | **0.116** (r > 0.3: 3개) | −0.188 | −0.319 |

- 무작위와의 차이: **H-optimus-0 +0.44**, Hibou-L +0.15, ResNet50 0. 타일 과제는 셋 다 세포 구성을 잘 담음 (R² 0.6–0.94), 유전자 과제에서 차이가 크게 벌어짐.
- 50명 · ridge 하나의 결과라 수치 자체는 흔들림. 다만 차이가 커서 순위는 믿을 만함.
- 속도: GPU 속도와 실제 추출의 읽기 속도(초당 약 640)로 계산. H-optimus-0도 슬라이드당 약 35초로 CellViT(약 22분)에 비해 작음.
- **결정: H-optimus-0 하나만 뽑음** (1536차원 float16, 510장 약 12 GB). 첫 run_batch 실행이 ResNet50으로 뽑은 150장은 서버에서 지움 (남겨 두면 "h5가 있으면 슬라이드 삭제" 규칙에 걸려 H-optimus-0 특징 없이 슬라이드가 지워질 위험. `extract.py`는 인코더가 다른 h5를 다시 뽑음).

**H-optimus-0으로 정한 이유 (요약)**
1. **목표 과제에서 유일하게 분명한 신호:** 단계 2는 H&E로 유전자 발현을 예측하는 것. 50명 축소판에서 무작위보다 r 중앙값이 +0.44 높고 R 목록 82개 유전자가 r > 0.3. ResNet50은 무작위와 같고(+0.01), Hibou-L은 약함(+0.15).
2. **세포 구성도 가장 잘 담음:** 타일 과제 R²가 대체로 가장 높고, 특히 염증 세포(0.82 vs ResNet50 0.70). 면역 세포는 재발 관련 유전자와 연결될 수 있음.
3. **속도 부담이 작음:** 510장 약 5시간, 슬라이드당 약 35초로 CellViT(장당 약 22분)에 비해 작아 전체 일정에 거의 영향 없음.
4. **라이선스 · 데이터 조건을 지킴:** Apache 2.0이라 특징 공개 가능, TCGA로 학습하지 않아 TCGA 평가 때 누수 없음, 동의만 하면 자동 승인.
5. **공개 벤치마크와 일치:** eva 9과제 평균 0.793으로 후보 중 가장 높고, Neidlinger et al. 2025 비교에서도 상위권.
6. **하나만 뽑는 이유:** ResNet50은 이 과제에서 신호가 없어 기준선으로서 의미가 작고, 두 개를 뽑으면 슬라이드를 오래 보관해야 함. 최종 성능은 B4에서 500명 교차검증으로 확인.
- **H-optimus-0 특징이 생기기 전에는 슬라이드를 지우지 않음.**

**슬라이드 QC 결과 (`slide_qc.py`, TCGA-BRCA 진단 슬라이드 1,133장, 2026-11-04)**

| 판정 | 전체 | 1,047명 | 515명 | 처리 |
|---|---|---|---|---|
| H&E | 1,062 | 978 | 483 | 사용 |
| H&E_faded (에오신 바램, 52장이 A8 기관) | 55 | 53 | 28 | **사용** (아래 시험 결과) |
| **IHC_suspect** | 1 | 1 | 1 | **제외** (TCGA-5L-AAT0) |
| check_thumb (GDC 썸네일 깨짐) / error (타일 API 400 · 500) | 5 / 11 | 4 / 11 | 1 / 3 | 원본 받아서 확인 |
| 펜 (초록 · 파랑 > 0.5%) | 12 | 11 | 9 | 펜 타일만 제외 |

- **바랜 H&E 시험:** TCGA-A8-A096 (H&E 색 비율 0.004, 가장 바램)과 A8-A0AD (0.46)를 CellViT로 돌림. 조직 타일당 핵 수 84 · 83개로 정상 H&E (88 · 90개)와 비슷하고, 타일 종양 비율 분포도 같은 모양 (0 근처와 1 근처에 몰림). QuPath 확대에서도 종양 · 기질 구분이 맞음 → 색 정규화 없이 사용.
- **빨간 표시 (`red_kind`):** 빨강 > 0.5%인 39장을 눈으로 나눈 결과에 맞춤 — 빨간 펜 선 16 (E9 기관, 펜 타일 제외), 지워진 흔적 15 (잉크 픽셀만 제외), 조직 자체가 빨강 7 (혈액 · 진한 에오신, 펜 아님), 썸네일 이상 1.
- **TCGA-5L-AAT0 제외:** GDC에 `Diagnostic Slide`(DX1)로 등록돼 있으나 실제로는 상피 표지자 IHC (청회색 + DAB 갈색, 에오신 없음, 2021-10 등록). 같은 기관 다른 환자(5L-AAT1)는 H&E → 이 슬라이드만의 업로드 오류로 봄. H&E 갈래(슬라이드 다운로드 · 타일 · 특징)에서 제외하고, RNA만 쓰는 분석에는 남김. H&E 갈래 515명 → **514명**.

**단계 2. "잘 예측된 유전자" 판정 (SEQUOIA와 동일)**
- 환자 단위 5-fold, 테스트 예측을 모두 합쳐 유전자마다 계산.
- ① Pearson r > 0, P < 0.05 ② 학습하지 않은 무작위 모델보다 r이 유의하게 큼 (Steiger, P < 0.05, BH < 0.2) ③ RMSE가 무작위 모델보다 작음.
- 순위: Pearson r 내림차순 (명시).
- 민감도 분석: 기관(TSS) 단위 fold, FDR < 0.05 + r ≥ 0.3.

**단계 3. 원격 재발 유전자**
- 사건 = DR, 시간 = TDR, 나머지는 마지막 추적일 또는 사망 시 중도 절단.
- 주 분석: 재발 기록 완전 환자. 민감도: 전체 1,413명, 내분비요법만 받은 환자.
- 학습은 HR+/HER2− 전체(병기 보정), 평가는 검사 권고 코호트.
- **Cox 회귀를 쓴 이유**
  - 재발 없이 추적이 끝난 환자(중도 절단)를 버리지 않고 "그때까지는 재발 없음"으로 씀 (1,251명 중 재발 272건뿐). t-검정 · 로지스틱은 추적 기간 차이를 무시함.
  - 발현을 높음/낮음으로 자르지 않고 연속값 그대로 넣음 (Kaplan–Meier · 로그순위 검정은 잘라야 함).
  - 림프절 · 크기 · 등급 · 나이를 같이 넣어 "임상 정보를 넘어서는 유전자"만 남길 수 있음.
  - 재발 위험이 시간에 따라 어떻게 변하는지 가정하지 않음 (와이블 등 모수 모형은 가정 필요).
  - 결과가 "1 SD당 위험비" 하나라 위험 / 보호 유전자가 바로 보이고, Oncotype · SEQUOIA 등 기존 서명 연구와 같은 방식이라 비교가 쉬움. 단계 5도 Cox라 이어짐.
  - 유전자 17,813개를 하나씩 돌려도 빠름.
  - 한계: 다른 원인 사망을 중도 절단으로 처리함 (재발 확률 자체를 보려면 Fine–Gray). 위험비가 추적 기간 전체 평균이라 초기 재발에 강한 증식 유전자는 효과가 섞임 (비례위험 가정 미확인).

**단계 5. 재발 수치**
- 유전자 발현은 코호트 안에서 유전자별 z-score로 맞춘 뒤 공식 적용 (마이크로어레이 → RNA-seq 이동 대비).
- 출력: Cox 선형 점수(0–100 변환) + METABRIC 기준 생존곡선으로 10년 원격 재발 확률.
- 위험 구간 기준은 METABRIC에서 정하고 고정, SCAN-B에서는 바꾸지 않음.

---

## 4. 다음 작업

두 갈래로 나눠 진행합니다. **A(RNA)는 슬라이드 없이 지금 바로** 할 수 있고, **B(H&E)는 슬라이드 저장 공간 문제를 먼저 풀어야** 합니다. C는 A와 B가 끝나야 합니다.

```
A (RNA, 지금)        A1 → A2 → A3 → A4 ──┐
                              └→ A5 ─────┤
B (H&E, 저장 공간)   B1 → B2 → B3 → B4 ──┼→ C1 → C2
```

### A. RNA 갈래 (METABRIC · SCAN-B)
| # | 작업 | 입력 → 출력 | 완료 기준 |
|---|---|---|---|
| **A1** | 발현 전처리 | METABRIC: 프로브 → 유전자 평균 (이미 log2). SCAN-B: Ensembl → 유전자 기호, log2(FPKM + 0.1). METABRIC ∩ SCAN-B ∩ TCGA 공통 유전자만 → `data/{metabric,scanb}/expr_hr_her2neg.parquet`. z-score는 분석 범위마다 `expression.zscore()`로 계산 | **완료**: 공통 **17,813개** 유전자 (Oncotype 21개 모두 포함). METABRIC 1,413명, SCAN-B 5,037명 (DRFi 없는 환자 포함). METABRIC 결측 13개 (Oncotype 외 11개 유전자, 1개씩) → 사용할 때 유전자 중앙값으로 채움 |
| **A2** | 코호트 표 고정 | scratchpad 집계 코드를 `cohort.py`로 옮김 (`standardize_metabric`, `standardize_scanb`, 기존 `assign_group` 재사용). Rueda DR/TDR, 재발 기록 완전, 내분비요법만 여부 → `data/{metabric,scanb}/cohort.csv` | **완료** (METABRIC 검사 권고 + 기록 완전 989명 / 원격 재발 179건, SCAN-B 검사 권고 3,095명 / 128건) |
| **A3** | 유사 RS 예후력 확인 | `pseudo_rs.py`를 METABRIC(마이크로어레이 log2 값) · SCAN-B(log2 FPKM)에 적용 → 10년 원격 재발 C-index, 순위 3구간 KM, 내분비요법만 하위군 | **완료** (`labels/scripts/pseudo_rs_prognosis.py`, 결과 `outputs/module1/a3_*`): 검사 권고 코호트에서 10년 원격 재발 C-index **METABRIC 0.625** (95% CI 0.59–0.66, 154건) · **SCAN-B 0.705** (0.66–0.75, 127건). 10년 재발률 low / intermediate / high = METABRIC 11.6 / 16.8 / 27.5%, SCAN-B 4.2 / 6.1 / 14.4%. 내분비요법만 하위군에서도 high가 갈림 (log-rank p < 0.001). 단, low와 intermediate의 차이는 내분비요법만 하위군에서 작음 (METABRIC 14.5 vs 15.6%, SCAN-B 4.9 vs 5.5%) |
| **A4** | **단계 3: R 목록** | METABRIC HR+/HER2− 기록 완전 1,251명 (원격 재발 272건). 유전자별 Cox (log 림프절 수 · 크기 · 등급 · 나이 보정, BH FDR) → FDR < 0.05 전부를 R 목록으로 → `outputs/module1/r_genes.csv` (`labels/scripts/r_genes.py`) | **완료**: **493개** (위험 332 / 보호 161). 보정 없이는 2,054개, 임상 보정으로 1,574개 빠짐. Oncotype 유전자는 증식 4개(AURKA · BIRC5 · CCNB1 · MKI67)만 남음. 처음 계획한 Cox-LASSO · elastic net 안정성 선택은 증식 유전자 묶음에서 대표만 남겨(선택 비율 0–24%) 교집합 목적에 맞지 않아 뺐음. SCAN-B 확인: 493개 Cox 점수 C-index 0.691, elastic net 29개 0.649, 무작위 29개 평균 0.612 |
| **A5** | 비교 기준선 준비 | 증식 점수 (Oncotype 증식 5개), NPI (METABRIC 열 / SCAN-B는 크기 · 림프절 · NHG로 계산), PAM50 ROR (SCAN-B 제공값), 무작위 유전자 세트, SEQUOIA 272개 서명 (Supplementary Data 5 확보 필요) | **부분 완료** (`scripts/baselines.py`, 결과 `data/{metabric,scanb}/baselines.csv`, `outputs/module1/a5_baselines.json`). 검사 권고 코호트 10년 원격 재발 C-index: METABRIC (989명) 유사 RS 0.625 · 증식 5개 **0.678** · SEQUOIA 0.580 · NPI 0.638 / SCAN-B (3,053명, 점수 모두 있는 환자) 유사 RS **0.702** · 증식 0.663 · SEQUOIA **0.702** · NPI 0.695 · PAM50 ROR 0.655. SEQUOIA는 272개 중 공통 유전자 259개, 계수 × (log2 발현 − 코호트 평균). METABRIC 값이 낮은 것은 TCGA RNA-seq 계수를 마이크로어레이에 그대로 쓴 영향으로 봄. 남은 것: CTS5 (공식 원문 확인), METABRIC ROR (제공값 없음), 무작위 유전자 세트 (C2에서 최종 서명과 같은 개수로) |

### B. H&E 갈래 (TCGA · CPTAC)
| # | 작업 | 입력 → 출력 | 완료 기준 · 주의 |
|---|---|---|---|
| **B1** | 슬라이드 받기 · QC | GPU 서버(`ssh gpu`, `data/sj/tcga/`)에 HR+/HER2− 514명 DX1 받기 (`slide_download.py`), 썸네일 QC (`slide_qc.py`) | **진행 중** (2026-11-05): 514장 중 510장 받음, 실패분 재시도 중. 1,047명 확장 · 나머지 532명 RNA는 필요할 때 |
| **B2** | **단계 1: 종양 영역** | 세포 분할 모델 **CellViT-SAM-H** (단계 1 세부 기준 참고) → 타일별 종양 핵 비율 → 기준값 이상 타일 목록 (`tumor_tiles.py`) | **기준 확정** (6장 규칙, 6-1 RNA 조각 맞춤). 514장 처리 진행 중 (GPU 1장 약 7일) |
| **B3** | 특징 추출 | `features/scripts/extract.py`: **H-optimus-0** (CLS, 1536차원), 128 µm 타일 (20배 256픽셀 해당)을 224로 줄여 입력. 조직 타일 전부 + `tumor_tile` 표시를 슬라이드별 `features/{슬라이드}.h5`에 저장 → 학습 때 종양 / 전체 타일, DX1만 / 모든 슬라이드를 고름. 묶음 실행: `features/scripts/run_batch.py` (받기 → CellViT → 종양 타일 → 특징) | 인코더는 시험 실행으로 H-optimus-0 확정 (2026-11-06, 단계 1 특징 인코더 참고). 시험 속도 초당 224타일 (RTX 5090), 510장 약 5시간. 종양 범위는 CellViT로 확정 (BCSS 비교 위 참고). **H-optimus-0 특징이 생긴 슬라이드만 지움** |
| **B4** | **단계 2: P 목록** | H&E 특징 → 유전자 발현 예측, SEQUOIA 기준으로 "잘 예측됨" 판정, Pearson r 순위 → `outputs/module1/p_genes.csv`. CPTAC (RNA 133명, 슬라이드 113 GB, ER/HER2는 Krug 2020 Table S1)로 외부 확인 | RNA 조각 맞춤 / 고정 기준 / 전체 타일 세 버전 결과 비교 (`tile_sets.select`) |

### C. 합치기
| # | 작업 | 완료 기준 |
|---|---|---|
| **C1** | **단계 4: 교집합** P ∩ R | 교집합 크기, 무작위로 기대되는 크기와 비교 (초기하 검정), 증식 · 면역 · 기질 유전자 비율 |
| **C2** | **단계 5: 재발 수치** METABRIC 교차검증으로 Cox 학습 → 가중치 · 위험 구간 고정 → SCAN-B 검사 권고 코호트에서 A3 · A5 기준선과 비교 | 10년 원격 재발 C-index, 위험 구간별 재발률, 임상 저위험 중 고위험군의 재발률 (검사 권고 근거), 보정 곡선 |

H&E 단계 평가(TCGA에서 H&E 예측 서명 → PFI, Howard · SEQUOIA 방식과 비교)는 C2 이후에 정합니다.

## 5. 확인하지 못한 것
- CPTAC-BRCA의 ER/HER2 상태와 슬라이드 종류 (FFPE / 냉동)
- 썸네일이 깨지거나(check_thumb 6) 타일 API가 실패한(error 5) 11장의 원본 상태 → `run_batch`가 CellViT를 돌리지 않고 목록만 출력
- 조직 접힘 · 연속 절편 판정 (접힘은 종양 핵도 진해서 단순 진하기 마스크는 종양을 지울 위험 → 보류)
- ~~SEQUOIA Supplementary Data 5 (272개 유전자 계수)의 공개 여부와 이용 조건~~ → 확보: 논문 MOESM4 엑셀 `Supp. Data 5` 시트 (`data/sequoia/raw/`), 논문 라이선스 CC BY-NC-ND 4.0이라 파일은 저장소에 올리지 않음
- SCAN-B 폐경 상태 (나이로 근사함)

---

## 6. B 갈래 처리 규칙 (종양 타일 · 특징 추출, 2026-11-05 확정)

| 단계 | 규칙 | 코드 |
|---|---|---|
| 대상 슬라이드 | HR+/HER2− manifest 환자의 첫 진단 슬라이드(DX1). `--all_slides`면 환자의 DX 전부 | `slide_download.py` |
| 슬라이드 제외 | QC `IHC_suspect` (TCGA-5L-AAT0) 받지 않음. `check_thumb` · `error` · QC 없음은 CellViT 보류 (원본 확인) | `slide_qc.py`, `run_batch.py` |
| 세포 분류 | CellViT-SAM-H, PanNuke 5종. 결과는 핵 중심 · 종류만 남김 (`cell_detection.json`) | `run_batch.py` |
| 타일 격자 | 128 µm 정사각형 (20배 256픽셀). 원본 픽셀 = round(128 / mpp) (0.252 mpp → 508, 0.2325 mpp → 551). 원점 (0, 0)부터 겹침 없이 | `tumor_tiles.py` |
| 펜 제외 | 초록 · 파랑 · 검정 펜 + 빨간 펜 선 · 흔적 (조직이 빨간 경우 제외)을 썸네일에서 5×5 최댓값 필터로 넓힌(사방 2픽셀) 마스크가 타일의 1% 넘게 겹치면 제외 | `slide_qc.py`, `tumor_tiles.py` |
| 조직 타일 | 핵 ≥ 20 & 펜 아님 | `tumor_tiles.py` |
| 종양 타일 (고정 기준) | 조직 타일 & Neoplastic 핵 / 전체 핵 ≥ **0.7** & Neoplastic 핵 ≥ **50** (2026-11-05 변경, 이전 0.5 · 개수 기준 없음). 0.7은 RNA 조각의 종양 핵 % (중앙값 80, 25% 70)에 맞춤, 50개는 세포가 듬성한 섬유질 위주 타일 제외 | `tumor_tiles.py`, `tile_sets.py` |
| 종양 타일 (RNA 조각 맞춤) | 종양 핵 ≥ 50인 조직 타일을 종양 핵 비율 높은 순으로 담다가, 담은 타일 전체의 종양 핵 비율(종양 핵 합 / 전체 핵 합)이 **환자의 RNA 조각 종양 핵 %**까지 내려오면 멈춤. 예측 때(새 환자는 값 없음)는 코호트 중앙값 80% | `rna_purity.py`, `tile_sets.py` |
| 보조 확인 | 큰 펜 원(볼록 껍질) 안쪽 종양 타일 비율, 썸네일 · 선택 결과 비교 JPG | `tumor_tiles.py` |
| 특징 | 조직 타일 전부를 원본에서 읽어 224×224로 줄여 **H-optimus-0** (ViT-g/14, CLS 출력) → 1536차원 float16. 정규화는 모델 카드 값 | `extract.py` |
| 저장 | 슬라이드마다 `features/{슬라이드}.h5`: features, coords [x, y], tumor_tile, tumor_ratio, n_cells, 속성 tile · mpp · encoder | `extract.py` |
| 학습 때 고르는 것 | RNA 조각 맞춤 / 고정 기준 / 조직 타일 전부 × DX1만 / 환자 슬라이드 전부 (환자 단위로 이어 붙인 뒤 고름). h5의 `tumor_tile`은 추출 당시 기준(0.5)이라 쓰지 않고 `tumor_ratio` · `n_cells`로 다시 판정 | `tile_sets.select` (B4) |
| 슬라이드 삭제 | H-optimus-0 특징(`features/{슬라이드}.h5`, encoder 속성 `h_optimus_0`)이 생긴 슬라이드만 | 서버에서 수동 |

색 정규화는 하지 않음 (바랜 H&E에서도 CellViT 핵 검출 · 종양 판정이 정상 H&E와 비슷했음).

### 6-1. 종양 타일 기준을 RNA 조각에 맞춘 이유 (2026-11-05)
- **RNA와 사진은 다른 조직 조각에서 나왔습니다.** TCGA는 냉동 조각에서 RNA를 뽑고, 그 조각의 위(TS) · 아래(BS) 단면을 병리의사가 판독해 종양 핵 % 등을 남겼습니다. 진단 슬라이드(DX, 파라핀)는 다른 조각이고 이 값이 없습니다 (1,130장 중 1장).
- **RNA 조각 판독값은 우리 515명 모두에게 있습니다** (TCGA-BRCA 전체 1,097 / 1,098명). TS만 304명, TS+BS 189명, BS만 19명, MS만 3명. 종양 핵 평균 76.6% (중앙값 80, 60% 미만 20명), 기질 22.2%, 정상 4.6%, 림프구 5.6%. 눈으로 어림한 5% 단위라 TS와 BS가 평균 8.3%p 다르고 189명 중 24명은 20%p 이상 다름 → ±10%p 정도로 맞추면 충분.
- **RNA는 세포에서 나오므로 면적이 아니라 핵 비율로 맞춥니다.** TCGA의 시료 기준(종양 핵 ≥ 60%)도 같은 지표.
- **DX 슬라이드 전체는 RNA 조각보다 종양이 훨씬 적습니다.** 처리된 15장에서 CellViT로 센 종양 핵 비율과 RNA 조각 값의 평균 차이:

| 타일 범위 | RNA 조각과의 평균 차이 (%p) |
|---|---|
| 조직 타일 전부 | 33.2 (DX 평균 약 45% vs RNA 조각 약 77%) |
| 고정 기준 (≥ 0.7 & 종양 핵 ≥ 50) | 11.0 (항상 80–92%로 RNA 조각보다 높음) |
| **RNA 조각 맞춤** | **5.6** (15명 중 9명은 목표값과 같음. 6명은 종양 핵 ≥ 50 후보를 다 담아도 목표보다 높아 후보 전부) |

- **타일마다 "비율 ≥ 목표"로 자르면 안 되는 이유:** 고른 타일 전체의 평균은 목표보다 높아짐. 평균을 직접 맞춰야 함.
- **모든 타일에 가중치만 주면 안 되는 이유:** 평균 숫자는 맞출 수 있지만, 멀리 떨어진 정상 조직 · 지방 경계가 섞임. RNA 조각의 비종양 부분은 주로 종양 사이 기질. 비율 높은 순으로 담으면 종양 중심에서 경계 쪽으로 넓혀 가는 macrodissection과 비슷해짐.
- 세 버전(RNA 조각 맞춤 · 고정 기준 · 전체 조직)을 B4에서 같이 돌려 비교. CellViT와 특징은 다시 뽑지 않음.
- CellViT가 센 핵 비율과 병리의사의 어림값은 같은 자가 아니므로, 위 비교는 방향을 보는 용도. 대량 처리 후 515명으로 다시 확인.

---

## 7. 실행 순서 (conda `team_project`, 프로젝트 폴더에서)

**A. RNA 갈래 (이 PC)**
| 순서 | 명령 | 결과 |
|---|---|---|
| 1 | `python -m module1_recurrence.data.scripts.external_fetch --expression` | METABRIC · SCAN-B 원본 |
| 2 | `python -m module1_recurrence.data.scripts.external_cohort` | `data/{metabric,scanb}/cohort.csv` |
| 3 | `python -m module1_recurrence.data.scripts.expression` | `data/{metabric,scanb}/expr_hr_her2neg.parquet` |
| 4 | `python -m module1_recurrence.labels.scripts.pseudo_rs_prognosis` | A3: `data/*/pseudo_rs.csv`, `outputs/module1/a3_*` |
| 5 | `python -m module1_recurrence.labels.scripts.r_genes` | A4: `outputs/module1/r_genes.csv` (R 목록), `r_gene_cox_all.csv` |
| 6 | `python -m module1_recurrence.scripts.baselines` | A5: `data/*/baselines.csv`, `outputs/module1/a5_baselines.json` |

SEQUOIA 파일은 `baselines.py` 맨 위 설명의 주소에서 받아 `data/sequoia/raw/`에 둔다.

**B. H&E 갈래**
| 순서 | 어디서 | 명령 | 결과 |
|---|---|---|---|
| 1 | 이 PC | `python -m module1_recurrence.data.scripts.gdc_fetch` → `...data.scripts.manifest` | `data/tcga/raw/`, `manifest.csv` |
| 2 | 이 PC | `python -m module1_recurrence.data.scripts.slide_qc` → `scp -r data/tcga/qc gpu:~/team_project_final/data/sj/tcga/` | 썸네일 · 펜 마스크 · `slide_qc.csv` |
| 3 | 서버 | `python -m module1_recurrence.data.scripts.slide_download --root data/sj/tcga --batch 600` | `slides/` |
| 4 | 서버 | `python -m module1_recurrence.features.scripts.run_batch --root data/sj/tcga --batch 0 --gpu 1` (반복) | `cellvit/`, `tiles/`, `features/` |
| 5 | 이 PC | `python -m module1_recurrence.data.scripts.rna_purity` | `data/tcga/rna_purity.csv` (환자별 RNA 조각 구성) |
| 6 | 서버 | 기준을 바꿨을 때 `python -m module1_recurrence.features.scripts.tumor_tiles --root data/sj/tcga --redo` | `tiles/` CSV · JPG 다시 계산 (특징은 그대로) |

---

## 8. 용어

| 용어 | 뜻 |
|---|---|
| **Cox 모델** | 사건(재발)이 언제 일어나는지를 다루는 회귀. 추적이 끝난 환자도 "그때까지는 재발 없음"으로 씀 |
| **위험비 (HR)** | 값이 1 단위 오를 때 재발 위험의 배수. 여기서는 유전자 발현 1 표준편차당 |
| **임상 보정** | 림프절 · 크기 · 등급 · 나이를 고정한 상태에서도 유전자가 추가로 연관되는지 봄 |
| **FDR** | 고른 목록 중 가짜 발견의 비율. FDR < 0.05면 493개 중 가짜는 평균 약 25개 이하 |
| **C-index** | 두 환자 중 먼저 재발한 환자에게 더 높은 점수를 줬을 확률. 0.5는 무작위, 0.6–0.7은 유전자 서명 수준 |
| **무작위 유전자 기준선** | 유방암에서는 아무 유전자나 골라도 증식과 연관돼 C-index 약 0.61이 나옴. 이 값을 넘어야 의미 있음 |
| **PQ · mPQ** | 핵 분할 지표. 탐지 F1 × 맞힌 핵의 평균 IoU (0–1). mPQ는 세포 종류별 PQ 평균 |
| **종양 타일** | 128 µm 타일 중 CellViT가 찾은 핵의 절반 이상이 Neoplastic인 타일 (6장) |
