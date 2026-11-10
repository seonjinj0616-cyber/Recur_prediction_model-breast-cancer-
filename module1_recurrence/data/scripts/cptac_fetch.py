"""CPTAC-BRCA (B4 외부 검증): GDC RNA-seq과 IDC 슬라이드(DICOM)를 받는다.

- RNA: GDC CPTAC-2 유방 원발 종양 STAR-Counts 133명 (TCGA와 같은 형식) → {root}/rnaseq/{patient}.tsv
- 슬라이드: IDC `cptac_brca` (TCIA CPTAC-BRCA의 DICOM 변환본, H&E 653장 · 198명) 중 RNA가 있는 환자의 슬라이드 전부.
  슬라이드 ID(ContainerIdentifier)가 RNA 샘플 ID와 같으면 RNA를 뽑은 그 조각의 슬라이드 → slides.csv의 rna_piece.
  슬라이드마다 폴더 하나 {root}/slides/{patient}_{slide_id}/*.dcm (CellViT · OpenSlide 4는 그 안의 .dcm 하나를 경로로 받는다).
이미 받은 것은 건너뛰므로 끊겨도 다시 실행하면 된다. `pip install idc-index`가 필요하다.
실행: python -m module1_recurrence.data.scripts.cptac_fetch --root data/sj/cptac [--rna_only] [--limit N]
결과: {root}/rna_files.csv, {root}/slides.csv, {root}/rnaseq/, {root}/slides/
"""
import argparse
import json
import shutil
import urllib.request
from pathlib import Path

import pandas as pd

from module1_recurrence.data.scripts.gdc_fetch import API, RAW, download


def rna_files():
    body = {"filters": {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": ["CPTAC-2"]}},
        {"op": "in", "content": {"field": "cases.primary_site", "value": ["Breast"]}},
        {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
        {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}}]},
        "fields": "file_id,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type", "size": 1000}
    req = urllib.request.Request(f"{API}/files", json.dumps(body).encode(), {"Content-Type": "application/json"})
    hits = json.load(urllib.request.urlopen(req))["data"]["hits"]
    rows = [{"patient_id": h["cases"][0]["submitter_id"], "rna_sample_id": h["cases"][0]["samples"][0]["submitter_id"],
             "sample_type": h["cases"][0]["samples"][0]["sample_type"], "rna_file_id": h["file_id"]} for h in hits]
    df = pd.DataFrame(rows).sort_values("patient_id")
    assert df.patient_id.is_unique, "환자당 RNA 파일이 하나가 아님"
    return df


def slide_list(client, rna):
    client.fetch_index("sm_index")
    sm = client.sql_query("""
        SELECT i.PatientID AS patient_id, s.ContainerIdentifier AS slide_id, i.SeriesInstanceUID AS series_uid,
               i.series_size_MB AS size_mb, s.ObjectiveLensPower AS objective, s.min_PixelSpacing_2sf AS spacing_mm,
               s.primaryAnatomicStructureModifier_CodeMeaning AS modifier
        FROM index i JOIN sm_index s USING (SeriesInstanceUID) WHERE i.collection_id = 'cptac_brca'""")
    sm = sm[sm.patient_id.isin(rna.patient_id)].copy()
    sm["slide"] = sm.patient_id + "_" + sm.slide_id
    sm["rna_piece"] = sm.slide_id.str[:30].isin(set(rna.rna_sample_id))
    return sm.sort_values("slide").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=RAW.parents[1] / "cptac", help="저장 폴더 (서버: data/sj/cptac)")
    ap.add_argument("--rna_only", action="store_true", help="RNA만 받기")
    ap.add_argument("--limit", type=int, help="슬라이드 앞에서 N장만 (테스트용)")
    args = ap.parse_args()

    (args.root / "rnaseq").mkdir(parents=True, exist_ok=True)
    rna = rna_files()
    rna.to_csv(args.root / "rna_files.csv", index=False)
    for pid, fid in zip(rna.patient_id, rna.rna_file_id):
        path = args.root / "rnaseq" / f"{pid}.tsv"
        if not path.exists():
            download(fid, path.with_suffix(".part"))
            path.with_suffix(".part").rename(path)
    print(f"RNA {len(rna)} patients → {args.root / 'rnaseq'}")
    if args.rna_only:
        return

    from idc_index import IDCClient
    client = IDCClient()
    sm = slide_list(client, rna)
    sm.to_csv(args.root / "slides.csv", index=False)
    print(f"slides {len(sm)} of {sm.patient_id.nunique()} RNA patients, {sm.size_mb.sum() / 1000:.1f} GB, "
          f"RNA piece {int(sm.rna_piece.sum())} slides / {sm[sm.rna_piece].patient_id.nunique()} patients")

    slide_dir = args.root / "slides"
    slide_dir.mkdir(exist_ok=True)
    todo = sm[[not (slide_dir / s).exists() for s in sm.slide]].head(args.limit)
    for i, r in enumerate(todo.itertuples(), 1):
        tmp = slide_dir / f".part_{r.slide}"
        shutil.rmtree(tmp, ignore_errors=True)  # 끊긴 다운로드는 처음부터
        client.download_from_selection(seriesInstanceUID=[r.series_uid], downloadDir=tmp, dirTemplate="", quiet=True)
        if not list(tmp.glob("*.dcm")):
            print(f"{r.slide} failed: no .dcm")
            continue
        tmp.rename(slide_dir / r.slide)  # 다 받은 슬라이드만 제 이름 폴더를 가진다
        print(f"{i}/{len(todo)} {r.slide} {r.size_mb:.0f} MB", flush=True)


if __name__ == "__main__":
    main()
