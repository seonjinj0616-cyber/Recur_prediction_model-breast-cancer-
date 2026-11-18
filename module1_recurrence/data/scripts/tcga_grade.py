"""TCGA HR+/HER2−: Nottingham 조직학적 등급 표 (CTS5 기준선용). TCGA 임상 파일에는 등급이 없다.

출처: Thennavan et al. Cell Genomics 2021 부록 mmc3.xlsx (TCGA BRCA Histologic features, 병리 전문가 판독, Heng et al. J Pathol 2017).
관 형성 · 핵 다형성 · 유사분열 점수(각 1~3)를 더해 3~5 → 1, 6~7 → 2, 8~9 → 3 등급.
판독이 없는 환자 (source 열에 표시)
- manual: TCGA-A2-A0EW는 관 형성 3 · 핵 다형성 3이고 유사분열만 빠져 2 또는 3 등급 → 3으로 지정 (2026-11-18 결정)
- Y-05: 세 점수가 모두 빈 환자는 virtual_hospital 규칙 Y-05로 채움. 아형별 등급 비율(R-02, 한국 KBCR,
  HR+/HER2− 0.268/0.528/0.204)에서 열쇠별 난수로 뽑는다 (virtual_hospital vh/common.rng, vh/tcga.synth_grade와 같은 계산).
  열쇠는 virtual_hospital 가상환자 번호(local/source_map.csv), 그곳에 없는 환자는 TCGA 바코드.
  실제 판독이 아니므로 민감도 분석(판독 있는 환자만)과 함께 본다.
실행: python -m module1_recurrence.data.scripts.tcga_grade
결과: data/tcga/grade.csv (patient_id, grade, source, rng_key)
"""
import hashlib
import io
import random
import zipfile
from pathlib import Path

import pandas as pd
import requests

from module1_recurrence.data.scripts.external_fetch import DATA
from module1_recurrence.data.scripts.gdc_fetch import RAW

MMC3 = RAW / "thennavan2021_mmc3.xlsx"
SUPP_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9028992/supplementaryFiles"
VH_SOURCE_MAP = Path(__file__).resolve().parents[4] / "virtual_hospital" / "local" / "source_map.csv"
OUT = DATA / "tcga" / "grade.csv"
COMPONENTS = ["Epithelial tubule formation", "Nuclear pleomorphism", "Mitosis"]
MANUAL = {"TCGA-A2-A0EW": 3}
R02_HR = [0.268, 0.528, 0.204]  # virtual_hospital rules/references.csv R-02 HR+/HER2-


def y05(key: str) -> int:
    """virtual_hospital Y-05: rng(key, "Y-05") = sha256("key|Y-05") 앞 16자리 seed."""
    x = random.Random(int(hashlib.sha256(f"{key}|Y-05".encode()).hexdigest()[:16], 16)).random()
    return 1 if x < R02_HR[0] else 2 if x < R02_HR[0] + R02_HR[1] else 3


def fetch():
    if MMC3.exists():
        return
    r = requests.get(SUPP_URL, timeout=300)
    r.raise_for_status()
    MMC3.write_bytes(zipfile.ZipFile(io.BytesIO(r.content)).read("mmc3.xlsx"))


def main():
    fetch()
    d = pd.read_excel(MMC3)
    d["patient_id"] = d["Sample CLID"].astype(str).str[:12]
    d = d.drop_duplicates("patient_id").set_index("patient_id")
    score = d[COMPONENTS].apply(lambda c: c.astype(str).str.extract(r"score = (\d)")[0].astype(float))
    total = score.sum(axis=1, min_count=3)
    read = pd.cut(total, [2, 5, 7, 9], labels=[1, 2, 3]).astype(float)

    pids = pd.read_csv(DATA / "tcga" / "manifest.csv").patient_id
    # 열쇠가 바뀌면 Y-05 값이 바뀌므로 번호표가 없으면 멈춘다 (서버에는 virtual_hospital이 없음 → grade.csv를 복사)
    vp = pd.read_csv(VH_SOURCE_MAP).set_index("source_id").vp_id
    rows = []
    for pid in pids:
        if pid in MANUAL:
            rows.append((pid, MANUAL[pid], "manual", ""))
        elif pd.notna(read.get(pid)):
            rows.append((pid, int(read[pid]), "thennavan2021", ""))
        else:
            key = vp.get(pid, pid)
            rows.append((pid, y05(key), "Y-05", key))
    out = pd.DataFrame(rows, columns=["patient_id", "grade", "source", "rng_key"])
    out.to_csv(OUT, index=False)
    print(out.source.value_counts().to_dict())
    print(out[out.source != "thennavan2021"].to_string(index=False))


if __name__ == "__main__":
    main()
