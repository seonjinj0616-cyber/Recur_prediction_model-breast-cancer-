"""TCGA-BRCA 진단 슬라이드 전체를 원본 없이 썸네일로 점검한다 (염색 종류 · 펜 표시).

GDC 슬라이드 뷰어용 타일 API(DeepZoom)에서 긴 변 1024픽셀 이하 썸네일만 받는다 (슬라이드당 약 1 MB).
- 염색: 조직 픽셀 중 H&E 색(보라~분홍) 비율이 낮으면 IHC 의심 → 눈으로 확인한 뒤 분석에서 뺀다.
- 펜: 초록 · 파랑 · 검정 펜 픽셀 비율과 마스크(png)를 남긴다 → 타일 단계에서 펜이 지나가는 타일을 뺀다.
실행: python -m module1_recurrence.data.scripts.slide_qc [--limit N] [--root DIR]
결과: {root}/qc/slide_qc.csv, {root}/qc/thumbs/*.png, {root}/qc/pen_masks/*.png, {root}/qc/ihc_suspects.png
"""
import argparse
import io
import json
import math
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from module1_recurrence.data.scripts.gdc_fetch import API, RAW

MAX_SIDE = 1024
RETRIES = 3
IHC_HE = 0.5  # 조직 중 H&E 색(보라~분홍) 비율이 이보다 낮으면 IHC 의심. 확인한 5장: IHC 0.06, H&E 0.63 이상


def _get(url):
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read()
        except OSError:
            if attempt == RETRIES - 1:
                raise
            time.sleep(5 * (attempt + 1))


def thumbnail(file_id):
    """DeepZoom 피라미드에서 긴 변이 MAX_SIDE 이하인 가장 큰 단계를 타일로 받아 이어 붙인다."""
    meta = json.loads(_get(f"{API}/tile/metadata/{file_id}"))
    w, h, ts, ov = int(meta["Width"]), int(meta["Height"]), int(meta["TileSize"]), int(meta["Overlap"])
    top = math.ceil(math.log2(max(w, h)))  # 원본 해상도 단계
    down = 2 ** max(0, math.ceil(math.log2(max(w, h) / MAX_SIDE)))
    level = top - int(math.log2(down))
    lw, lh = math.ceil(w / down), math.ceil(h / down)
    img = Image.new("RGB", (lw, lh), "white")
    for ty in range(math.ceil(lh / ts)):
        for tx in range(math.ceil(lw / ts)):
            tile = Image.open(io.BytesIO(_get(f"{API}/tile/{file_id}?level={level}&x={tx}&y={ty}"))).convert("RGB")
            # 첫 행 · 열이 아니면 타일 앞쪽에 겹침(overlap) 픽셀이 붙어 있다
            img.paste(tile, (tx * ts - (ov if tx else 0), ty * ts - (ov if ty else 0)))
    return img, down


def measure(img):
    """HSV 색 범위로 조직 · 분홍 · 갈색 · 펜 픽셀을 센다. 비율은 조직(펜 제외) 픽셀 기준."""
    hsv = np.asarray(img.convert("HSV")).astype(int)
    hue, sat, val = hsv[..., 0] * 360 // 255, hsv[..., 1], hsv[..., 2]
    green = (hue >= 90) & (hue <= 200) & (sat > 60)
    blue = (hue >= 200) & (hue <= 260) & (sat > 110) & (val < 200)
    dark = val < 60
    pen = green | blue | dark
    tissue = (sat > 20) & (val < 235) & ~pen
    n = max(int(tissue.sum()), 1)
    # H&E 조직은 보라~분홍(270–360°, 0–15°), IHC는 헤마톡실린만 있어 청회색(200–270°) + 갈색 DAB
    he = tissue & ((hue >= 270) | (hue <= 15))
    hblue = tissue & (hue >= 200) & (hue < 270)
    brown = tissue & (hue > 15) & (hue <= 45) & (sat > 50)
    stats = {
        "tissue_frac": tissue.mean(),
        "he_frac": he.sum() / n,
        "hblue_frac": hblue.sum() / n,
        "brown_frac": brown.sum() / n,
        "green_pen_frac": green.mean(),
        "blue_pen_frac": blue.mean(),
        "dark_frac": dark.mean(),
    }
    return stats, pen


def run_one(row, out):
    name = row.file_name.split(".")[0]  # 예: TCGA-A1-A0SB-01Z-00-DX1
    thumb_path = out / "thumbs" / f"{name}.png"
    try:
        if thumb_path.exists():
            img, down = Image.open(thumb_path).convert("RGB"), None
        else:
            img, down = thumbnail(row.file_id)
            img.save(thumb_path)
        stats, pen = measure(img)
        Image.fromarray(pen.astype(np.uint8) * 255).save(out / "pen_masks" / f"{name}.png")
        return {"file_id": row.file_id, "slide": name, "patient_id": row.patient_id,
                "thumb_w": img.width, "thumb_h": img.height, "downsample": down, **stats}
    except OSError as e:
        print(f"{name} failed: {e}")
        return {"file_id": row.file_id, "slide": name, "patient_id": row.patient_id, "error": str(e)}


def contact_sheet(paths, dest, cell=256, cols=8):
    """IHC 의심 슬라이드 썸네일을 한 장에 모은다 (눈으로 확인용)."""
    rows = math.ceil(len(paths) / cols)
    sheet = Image.new("RGB", (cols * cell, max(rows, 1) * cell), "white")
    for i, p in enumerate(paths):
        im = Image.open(p)
        im.thumbnail((cell - 4, cell - 4))
        sheet.paste(im, ((i % cols) * cell + 2, (i // cols) * cell + 2))
    sheet.save(dest)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=RAW.parent, help="raw/dx_slides.tsv가 있는 폴더")
    ap.add_argument("--limit", type=int, help="앞에서 N장만 (테스트용)")
    args = ap.parse_args()

    out = args.root / "qc"
    for d in ("thumbs", "pen_masks"):
        (out / d).mkdir(parents=True, exist_ok=True)

    # 환자 단위로 거르지 않고 TCGA-BRCA 진단 슬라이드 전부 (1,047명 코호트를 포함)
    dx = pd.read_csv(args.root / "raw" / "dx_slides.tsv", sep="\t").rename(
        columns={"cases.0.submitter_id": "patient_id"}).sort_values("file_name").head(args.limit)

    with ThreadPoolExecutor(max_workers=8) as ex:
        rows = list(ex.map(lambda r: run_one(r, out), dx.itertuples()))

    qc = pd.DataFrame(rows)
    qc["ihc_suspect"] = qc.he_frac < IHC_HE
    qc.to_csv(out / "slide_qc.csv", index=False)

    suspects = qc[qc.ihc_suspect].sort_values("he_frac")
    contact_sheet([out / "thumbs" / f"{s}.png" for s in suspects.slide], out / "ihc_suspects.png")
    print(f"slides {len(qc)} | failed {int(qc.get('error', pd.Series()).notna().sum())} "
          f"| IHC suspects {len(suspects)} | pen (green+blue > 0.5%) "
          f"{int(((qc.green_pen_frac + qc.blue_pen_frac) > 0.005).sum())}")
    print(suspects[["slide", "he_frac", "hblue_frac", "brown_frac"]].to_string(index=False))


if __name__ == "__main__":
    main()
