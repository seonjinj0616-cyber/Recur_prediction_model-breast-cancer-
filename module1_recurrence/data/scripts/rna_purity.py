"""RNA를 뽑은 냉동 조직 조각의 구성(종양 핵 % 등)을 GDC에서 받는다 (B2, 환자별 맞춤 타일 선택용).

TCGA는 냉동 조직 조각에서 RNA를 뽑고, 그 조각의 위(TS) · 아래(BS) · 가운데(MS) 단면을 병리의사가 판독했다.
진단 슬라이드(DX, 파라핀)와는 다른 조직 조각이라 DX에는 이 값이 없다.
manifest의 RNA sample ID(예: TCGA-A1-A0SB-01A)와 같은 sample의 TS · BS · MS 판독값을 평균한다.
값은 병리의사가 눈으로 어림한 5% 단위이고, TS와 BS가 평균 8%p 다르다 (2026-11-05, 189명).

실행: python -m module1_recurrence.data.scripts.rna_purity
결과: data/tcga/rna_purity.csv (환자마다 한 줄, 값은 % 단위)
"""
import json
import urllib.request

import pandas as pd

from module1_recurrence.data.scripts.gdc_fetch import API, RAW, _eq

FIELDS = {"percent_tumor_nuclei": "tumor_nuclei", "percent_tumor_cells": "tumor_cells",
          "percent_stromal_cells": "stromal_cells", "percent_normal_cells": "normal_cells",
          "percent_lymphocyte_infiltration": "lymphocytes", "percent_necrosis": "necrosis"}


def fetch_sections():
    """TCGA-BRCA 모든 환자의 냉동 절편(TS · BS · MS) 판독값 → 절편마다 한 줄."""
    fields = ["submitter_id", "samples.submitter_id"] + \
             [f"samples.portions.slides.{f}" for f in ["submitter_id", *FIELDS]]
    body = {"filters": _eq("project.project_id", "TCGA-BRCA"), "fields": ",".join(fields),
            "format": "JSON", "size": "5000"}
    req = urllib.request.Request(f"{API}/cases", data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req) as r:
        hits = json.load(r)["data"]["hits"]
    rows = [{"patient_id": c["submitter_id"], "sample_id": s["submitter_id"], "section": sl["submitter_id"],
             **{v: sl.get(k) for k, v in FIELDS.items()}}
            for c in hits for s in c.get("samples", []) for p in s.get("portions") or []
            for sl in p.get("slides") or []]
    df = pd.DataFrame(rows)
    df["kind"] = df.section.str.extract(r"-(TS|BS|MS)\w*$")[0]  # DX(진단 슬라이드)는 빠짐
    return df[df.kind.notna() & df.tumor_nuclei.notna()]


def main():
    m = pd.read_csv(RAW.parent / "manifest.csv")
    sec = fetch_sections()
    sec = sec[sec.sample_id.isin(m.rna_sample_id)]
    out = sec.groupby(["patient_id", "sample_id"]).agg(
        sections=("kind", lambda k: "+".join(sorted(set(k)))),**{v: (v, "mean") for v in FIELDS.values()}).reset_index()
    out = m[["patient_id"]].merge(out, on="patient_id", how="left")
    out.to_csv(RAW.parent / "rna_purity.csv", index=False, float_format="%.1f")
    print(f"{out.tumor_nuclei.notna().sum()}/{len(m)} patients with RNA-piece tumor nuclei %")
    print(out.tumor_nuclei.describe().round(1).to_string())


if __name__ == "__main__":
    main()
