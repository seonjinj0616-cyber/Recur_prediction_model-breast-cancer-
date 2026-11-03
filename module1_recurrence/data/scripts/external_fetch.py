"""METABRIC · SCAN-B 원천 데이터를 받는다 (모두 공개, 가입 불필요).

METABRIC: cBioPortal brca_metabric 임상 · mRNA 샘플 목록, Rueda et al. Nature 2019 부록 Table S5/S6 (원격 재발 DR, TDR)
SCAN-B: Staaf et al. 2022 Mendeley Data (doi 10.17632/yzxtxn4nmd.3, CC BY 4.0) 임상 표 · 유전자 주석, GEO GSE96058 series matrix
발현 행렬은 커서(METABRIC 약 0.7 GB, SCAN-B 약 2.4 GB) --expression 을 줄 때만 받는다.
이미 있는 파일은 건너뛴다. --check 는 받지 않고 주소가 살아 있는지와 크기만 확인한다.
실행: python -m module1_recurrence.data.scripts.external_fetch [--expression] [--check]
결과: data/metabric/raw/, data/scanb/raw/, data/metabric/hr_her2neg.csv
"""
import argparse
import json
import urllib.request
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parents[3] / "data"
METABRIC = DATA / "metabric" / "raw"
SCANB = DATA / "scanb" / "raw"

CBIO = "https://www.cbioportal.org/api"
RUEDA = "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-019-1007-8/MediaObjects/41586_2019_1007_MOESM{}_ESM.txt"
MENDELEY = "https://data.mendeley.com/public-files/datasets/yzxtxn4nmd/files/{}/file_downloaded"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/matrix/{}"

FILES = {
    METABRIC / "clinical_patient.json": f"{CBIO}/studies/brca_metabric/clinical-data?clinicalDataType=PATIENT",
    METABRIC / "clinical_sample.json": f"{CBIO}/studies/brca_metabric/clinical-data?clinicalDataType=SAMPLE",
    METABRIC / "mrna_sample_ids.json": f"{CBIO}/sample-lists/brca_metabric_mrna/sample-ids",
    METABRIC / "rueda2019_TableS5.txt": RUEDA.format(7),  # 3,240명 전체
    METABRIC / "rueda2019_TableS6.txt": RUEDA.format(8),  # 발현 있는 1,980명
    SCANB / "scanb_clinical_2023-01-13.xlsx": MENDELEY.format("33d08b30-4685-4814-ab87-606a20c3092b"),
    SCANB / "Gene.ID.ann.txt": MENDELEY.format("4d4f41c8-efbd-4d3e-abb0-a210701dd634"),
    SCANB / "GSE96058-GPL11154_series_matrix.txt.gz": GEO.format("GSE96058-GPL11154_series_matrix.txt.gz"),
    SCANB / "GSE96058-GPL18573_series_matrix.txt.gz": GEO.format("GSE96058-GPL18573_series_matrix.txt.gz"),
}

EXPRESSION = {
    METABRIC / "data_mrna_illumina_microarray.txt":
        "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_metabric/data_mrna_illumina_microarray.txt",
    # 라이브러리 프로토콜 보정 FPKM (TruSeq 기준으로 맞춘 값)
    SCANB / "SCANB.9206.genematrix_noNeg.txt": MENDELEY.format("29d59830-01a6-49aa-82aa-8312ae669863"),
}

HEADERS = {"User-Agent": "Mozilla/5.0"}  # Mendeley가 기본 urllib 요청을 막는다


def fetch(url, path):
    tmp = path.with_suffix(path.suffix + ".part")
    with urllib.request.urlopen(urllib.request.Request(url, headers=HEADERS)) as r, open(tmp, "wb") as f:
        while chunk := r.read(1 << 20):
            f.write(chunk)
    tmp.rename(path)  # 다 받은 파일만 최종 이름을 가진다


def remote_size(url):
    req = urllib.request.Request(url, headers=HEADERS, method="HEAD")
    with urllib.request.urlopen(req) as r:
        return r.headers.get("Content-Length")


def build_hr_her2neg():
    """cBioPortal 임상 JSON → mRNA가 있는 HR+/HER2− 환자 표 (환자 1행)."""
    def wide(name, key):
        rows = json.loads((METABRIC / name).read_text(encoding="utf-8"))
        return pd.DataFrame(rows).pivot(index=key, columns="clinicalAttributeId", values="value")

    patient = wide("clinical_patient.json", "patientId")
    sample = wide("clinical_sample.json", "sampleId")  # METABRIC은 환자당 샘플 1개, sampleId = patientId
    df = patient.join(sample, how="inner")
    df.index.name = "patientId"

    mrna = set(json.loads((METABRIC / "mrna_sample_ids.json").read_text(encoding="utf-8")))
    hr = (df.ER_STATUS == "Positive") | (df.PR_STATUS == "Positive")
    df = df[hr & (df.HER2_STATUS == "Negative") & df.index.isin(mrna)]
    df.to_csv(DATA / "metabric" / "hr_her2neg.csv")
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--expression", action="store_true", help="큰 발현 행렬도 받기")
    ap.add_argument("--check", action="store_true", help="받지 않고 주소와 크기만 확인")
    args = ap.parse_args()

    files = FILES | (EXPRESSION if args.expression else {})
    for path, url in files.items():
        if args.check:
            local = f"{path.stat().st_size:,}" if path.exists() else "없음"
            print(f"{path.name}: 원격 {remote_size(url)} / 로컬 {local}")
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.exists():
            fetch(url, path)
        print(f"{path.name}: {path.stat().st_size:,} bytes")

    if not args.check:
        print("METABRIC HR+/HER2− with mRNA:", len(build_hr_her2neg()))


if __name__ == "__main__":
    main()
