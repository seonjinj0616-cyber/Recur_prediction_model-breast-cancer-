"""manifest.csv의 rna_file_id로 STAR-Counts RNA-seq 파일을 받는다.

먼저 data.manifest를 실행한다. 이미 받은 파일은 건너뛰므로 중간에 끊겨도 다시 실행하면 된다.
실행: python -m module1_recurrence.data.rnaseq_download [--limit N]
결과: data/tcga/rnaseq/{patient_id}.tsv
"""
import argparse
import time
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from module1_recurrence.data.gdc_fetch import RAW, download

MANIFEST = RAW.parent / "manifest.csv"
OUT_DIR = RAW.parent / "rnaseq"
RETRIES = 3


def fetch(patient_id, file_id):
    path = OUT_DIR / f"{patient_id}.tsv"
    if path.exists():
        return "skip"
    tmp = path.with_suffix(".part")
    for attempt in range(RETRIES):
        try:
            download(file_id, tmp)
            tmp.rename(path)  # 다 받은 파일만 .tsv 이름을 가진다
            return "ok"
        except OSError as e:  # GDC가 가끔 연결을 끊는다
            print(f"{patient_id} retry {attempt + 1}: {e}")
            time.sleep(5 * (attempt + 1))
    return "failed"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, help="앞에서 N명만 받기 (테스트용)")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    m = pd.read_csv(MANIFEST).dropna(subset=["rna_file_id"]).head(args.limit)

    with ThreadPoolExecutor(max_workers=4) as ex:
        results = list(ex.map(fetch, m.patient_id, m.rna_file_id))

    print({r: results.count(r) for r in set(results)}, "| files", len(list(OUT_DIR.glob("*.tsv"))))


if __name__ == "__main__":
    main()
