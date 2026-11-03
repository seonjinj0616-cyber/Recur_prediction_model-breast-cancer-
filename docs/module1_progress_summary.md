# 모듈 1 진행 정리 (RNA 갈래)

작성일: 2026-11-03 · 마지막 커밋 `ebb795c`
자세한 계획과 표는 `docs/module1_cohorts_plan.md`, 데이터 조사는 `docs/recurrence_validation_data.md`

---

## 1. 목표와 방향

**무엇을 만드나:** HR+/HER2− 조기 유방암에서 H&E 슬라이드로 원격 재발 위험을 알려 주는 모듈.

**위치:** "Oncotype DX 대체"가 아니라, **임상 정보로는 놓치기 쉬운 고위험 환자에게 Oncotype 검사를 권고하는 근거를 주는 모델.**

**차별화 5단계**
| 단계 | 내용 | 데이터 | 담당 |
|---|---|---|---|
| 1 | 종양이 밀집한 영역만 사용 | TCGA 슬라이드 | GPU 세션 |
| 2 | H&E로 잘 예측되는 유전자 (**P 목록**, SEQUOIA와 같은 기준) | TCGA H&E + RNA | GPU 세션 |
| 3 | 원격 재발과 연관된 유전자 (**R 목록**) | METABRIC | **완료 (이 세션)** |
| 4 | P ∩ R 교집합 | — | P 목록 대기 |
| 5 | 교집합 유전자로 재발 수치 계산, SCAN-B에서 기준선과 비교 | METABRIC 학습 → SCAN-B 평가 | 대기 |

핵심 아이디어: H&E로 보이는 유전자(P)와 재발에 중요한 유전자(R)를 **서로 다른 코호트에서 따로** 구한 뒤 겹치는 것만 쓴다.

---

## 2. 정한 원칙

| 항목 | 결정 |
|---|---|
| 검사 대상 (NCCN) | HR+/HER2−, N0이면서 종양 > 0.5 cm, 또는 N1(1–3개)이면서 폐경 후. pN1mi(미세전이)는 N0로 봄 |
| 종점 | **원격 재발** (국소 재발 제외). 평가는 10년에서 자름 |
| 학습 · 평가 | 학습은 HR+/HER2− 전체, 주 평가는 검사 권고 코호트 |
| 외부 검증 | 최종 비교는 **SCAN-B** (RNA-seq, 독립 코호트) |
| 데이터 공유 | GitHub에는 코드와 문서만. 데이터는 라이선스 · 용량 때문에 올리지 않음 |
| 환경 | conda `team_project` 환경에서 실행 (`conda activate team_project`) |

---

## 3. 데이터

| 데이터 | 내용 | 환자 수 | 원격 재발 |
|---|---|---|---|
| TCGA-BRCA | H&E + RNA-seq. 재발은 PFI뿐 (국소 + 원격 섞임, 추적 약 2년) | HR+/HER2− 515명 중 검사 권고 379명 | — |
| METABRIC | 마이크로어레이 + 장기 추적 (Rueda 2019 원격 재발 열) | HR+/HER2− 1,413명, 기록 완전 1,251명, 검사 권고 989명 | 272건 (검사 권고 154건, 10년 내) |
| SCAN-B | RNA-seq + 원격 재발(DRFi), CC BY 4.0 | HR+/HER2− 추적 코호트, 검사 권고 3,095명 | 128건 |
| SEQUOIA 서명 | 논문 부록 272개 유전자 계수 (CC BY-NC-ND 4.0) | — | — |

세 데이터 공통 유전자: **17,813개** (Oncotype 21개 모두 포함)

---

## 4. 한 일과 결과

| # | 작업 | 결과 | 커밋 |
|---|---|---|---|
| 1 | 유사 RS 구현 (Oncotype 공식을 RNA-seq에 적용) | Howard 연구용 RS와 순위 상관 0.86. 테스트 6개 통과 | `9447860` |
| 2 | 코호트 표 (A2) | METABRIC · SCAN-B 검사 권고 코호트 확정 | `9447860` |
| 3 | 발현 전처리 (A1) | 공통 유전자 17,813개로 정리 | `9d514a3` |
| 4 | 폴더 정리 | 폴더마다 `scripts/`에 .py 모음 | `9d514a3` |
| 5 | **유사 RS 예후력 (A3)** | 외부 두 코호트에서 원격 재발을 나눔 (아래) | `50e5e5d` |
| 6 | **R 목록 (A4)** | 493개 (아래) | `35871c2` |
| 7 | **비교 기준선 (A5)** | 같은 환자에서 기준선 5개 비교 (아래) | `ebb795c` |

### A3. 유사 RS가 원격 재발을 나누는가 → 나눈다
검사 권고 코호트, 10년 원격 재발

| | METABRIC (989명) | SCAN-B (3,095명) |
|---|---|---|
| C-index | 0.625 | 0.705 |
| 10년 재발률 low / intermediate / high | 11.6 / 16.8 / **27.5%** | 4.2 / 6.1 / **14.4%** |

- high 구간은 확실히 갈림 (내분비요법만 받은 환자에서도).
- low와 intermediate는 거의 차이가 없음 → 이 라벨의 정보는 대부분 "high인가 아닌가".

### A4. R 목록 = 493개
- 방법: 유전자 17,813개 각각에 Cox 모델 (림프절 수 · 크기 · 등급 · 나이 보정) → **FDR < 0.05** 전부.
- 결과: 위험 332개 / 보호 161개. 가장 강한 것은 증식 유전자 (AURKA, KIF20A, CCNB2, UBE2C, MKI67).
- 임상 보정 없이는 2,054개 → 보정으로 1,574개 빠짐 (림프절 · 등급과 같은 정보라서).
- Oncotype 유전자 중 남은 것: AURKA, BIRC5, CCNB1, MKI67 (증식 4개).
- **바꾼 결정:** 처음 계획한 LASSO · elastic net 선택은 비슷한 증식 유전자 묶음에서 대표 하나만 남겨서 뺐음. 교집합에는 "연관된 유전자 전부"가 필요하기 때문.
  - SCAN-B 확인: 493개 0.691 > elastic net 29개 0.649 > 무작위 29개 평균 0.612

### A5. 비교 기준선 (같은 환자, 10년 원격 재발 C-index)
| 점수 | METABRIC 검사 권고 | SCAN-B 검사 권고 |
|---|---|---|
| 유사 RS | 0.625 | **0.702** |
| SEQUOIA 272개 서명 | 0.580 | **0.702** |
| 증식 유전자 5개 | **0.678** | 0.663 |
| NPI (임상 지표) | 0.638 | 0.695 |
| PAM50 ROR | — | 0.655 |

- **우리 최종 모델의 목표:** SCAN-B 검사 권고 코호트에서 **0.70을 넘거나**, 기존 점수와 결합했을 때 추가 정보를 보여 주기.
- SEQUOIA가 METABRIC에서 낮은 것은 RNA-seq 계수를 마이크로어레이에 그대로 쓴 영향으로 봄 → 최종 비교는 SCAN-B 기준.

---

## 5. 용어 정리

| 용어 | 뜻 |
|---|---|
| **Cox 모델** | 사건(재발)이 언제 일어나는지를 다루는 회귀. 추적이 끝난 환자도 "그때까지는 재발 없음"으로 씀 |
| **위험비 (HR)** | 값이 1 단위 오를 때 재발 위험의 배수. 1.5면 1.5배, 0.7이면 30% 낮음. 여기서는 유전자 발현 1 표준편차당 |
| **임상 보정** | 림프절 · 크기 · 등급 · 나이를 고정한 상태에서도 유전자가 추가로 연관되는지 봄 |
| **FDR** | 고른 목록 중 가짜 발견의 비율. FDR < 0.05면 493개 중 가짜는 평균 약 25개(493 × 0.05) 이하로 예상. 실제 기준선은 p ≈ 0.0014 |
| **다변수 선택 (LASSO · elastic net)** | 여러 유전자를 한 모델에 넣고 꼭 필요한 최소 조합만 남김. 비슷한 유전자 중 대표만 남기므로 R 목록에는 맞지 않음. 5단계에서 과적합 방지용으로 씀 |
| **C-index** | 두 환자 중 먼저 재발한 환자에게 더 높은 점수를 줬을 확률. 0.5는 무작위, 0.6–0.7은 유전자 서명 수준, 0.7–0.8은 좋음 |
| **무작위 유전자 기준선** | 유방암에서는 아무 유전자나 골라도 증식과 연관돼 C-index 약 0.61이 나옴. 그래서 0.5가 아니라 이 값을 넘어야 의미 있음 |

---

## 6. 코드와 결과 파일

**실행 (conda `team_project` 환경, 프로젝트 폴더에서)**
| 순서 | 명령 | 결과 |
|---|---|---|
| 1 | `python -m module1_recurrence.data.scripts.external_fetch --expression` | METABRIC · SCAN-B 원본 받기 |
| 2 | `python -m module1_recurrence.data.scripts.external_cohort` | `data/{metabric,scanb}/cohort.csv` |
| 3 | `python -m module1_recurrence.data.scripts.expression` | `data/{metabric,scanb}/expr_hr_her2neg.parquet` |
| 4 | `python -m module1_recurrence.labels.scripts.pseudo_rs_prognosis` | A3: `data/*/pseudo_rs.csv`, `outputs/module1/a3_*` |
| 5 | `python -m module1_recurrence.labels.scripts.r_genes` | A4: `outputs/module1/r_genes.csv` (R 목록), `r_gene_cox_all.csv` |
| 6 | `python -m module1_recurrence.scripts.baselines` | A5: `data/*/baselines.csv`, `outputs/module1/a5_baselines.json` |

SEQUOIA 파일은 `baselines.py` 맨 위 설명의 주소에서 받아 `data/sequoia/raw/`에 둔다.

---

## 7. 남은 일

| 항목 | 상태 · 필요한 것 |
|---|---|
| **4단계 P ∩ R** | GPU 세션의 P 목록 파일 필요 |
| **5단계 재발 수치 (C2)** | 교집합 유전자로 METABRIC에서 Cox 학습 → SCAN-B에서 기준선과 비교, 무작위 유전자 세트 기준선 포함, 보정 곡선 확인 |
| A5 선택 항목 | CTS5 (공식 원문 확인), METABRIC PAM50 ROR (직접 계산) |
| 확인할 것 | METABRIC에서 증식 5개가 유사 RS보다 나은 이유, MYBL2가 METABRIC에서 증식 신호를 안 보이는 이유 (프로브 문제 의심) |
| H&E 갈래 저장 공간 | TCGA 슬라이드 약 1.09 TB 필요 (GPU 세션 쪽에서 처리) |
