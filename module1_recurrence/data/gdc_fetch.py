"""모듈 1 원천 데이터를 NCI GDC에서 받는다 (open access, 가입 불필요).

실행: python -m module1_recurrence.data.gdc_fetch
결과: data/tcga/raw/ 아래에 임상, TCGA-CDR, 슬라이드 목록, RNA-seq 목록
"""
import json
import urllib.request
from pathlib import Path

API = "https://api.gdc.cancer.gov"
RAW = Path(__file__).resolve().parents[2] / "data" / "tcga" / "raw"

# TCGA-CDR (Liu et al., Cell 2018) 부록. GDC PanCanAtlas 출판 페이지의 파일
CDR_FILE_ID = "1b5f413e-a8d1-4d10-92eb-7c4ae739ed81"
CLINICAL_FILE_NAME = "nationwidechildrens.org_clinical_patient_brca.txt"


def _eq(field, value):
    return {"op": "=", "content": {"field": field, "value": value}}


def query_files(filters, fields):
    """GDC /files 검색 결과를 TSV 문자열로 받는다."""
    body = {
        "filters": {"op": "and", "content": [_eq("cases.project.project_id", "TCGA-BRCA"), *filters]},
        "fields": ",".join(fields),
        "format": "TSV",
        "size": "10000",
    }
    req = urllib.request.Request(f"{API}/files", data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req) as r:
        return r.read().decode()


def download(file_id, path):
    with urllib.request.urlopen(f"{API}/data/{file_id}") as r:
        path.write_bytes(r.read())


def main():
    RAW.mkdir(parents=True, exist_ok=True)

    (RAW / "dx_slides.tsv").write_text(query_files(
        [_eq("data_type", "Slide Image"), _eq("experimental_strategy", "Diagnostic Slide")],
        ["file_id", "file_name", "file_size", "cases.submitter_id"]), encoding="utf-8")

    (RAW / "rnaseq.tsv").write_text(query_files(
        [_eq("data_type", "Gene Expression Quantification"),
         _eq("analysis.workflow_type", "STAR - Counts"),
         _eq("cases.samples.sample_type", "Primary Tumor")],
        ["file_id", "file_name", "file_size", "cases.submitter_id", "cases.samples.submitter_id"]),
        encoding="utf-8")

    # 임상 파일 UUID는 GDC 릴리스마다 바뀔 수 있어서 파일 이름으로 찾는다
    clinical = query_files([_eq("file_name", CLINICAL_FILE_NAME)], ["file_id"])
    clinical_id = clinical.splitlines()[1].split("\t")[0]
    download(clinical_id, RAW / CLINICAL_FILE_NAME)
    download(CDR_FILE_ID, RAW / "TCGA-CDR-SupplementalTableS1.xlsx")

    for f in sorted(RAW.iterdir()):
        print(f"{f.name}: {f.stat().st_size:,} bytes")


if __name__ == "__main__":
    main()
