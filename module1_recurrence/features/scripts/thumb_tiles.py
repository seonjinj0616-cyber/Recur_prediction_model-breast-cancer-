"""CellViT 없이 QC 썸네일로 조직 타일만 고른다 (C2 방법 4 · 5: 추가 학습 환자 · 외부 슬라이드용).

tumor_tiles와 같은 128 µm 격자를 쓰고, 타일 안 조직 픽셀(slide_qc.measure의 조직 기준: 채도 > 20, 밝기 < 235, 펜 제외) 비율이
MIN_TISSUE 이상이고 펜 비율이 PEN_MAX 이하인 타일을 조직 타일로 둔다. 핵 수 · 종양 비율은 없으므로 비워 둔다
(n_cells 0, tumor_ratio NaN, tumor_tile False) → extract.py가 그대로 읽고, 분석은 전체 조직(all) 방식만 쓸 수 있다.
썸네일이 {root}/qc/thumbs/에 없으면 원본에서 만든다 (긴 변 1024 이하, 2의 거듭제곱 축소).
실행: python -m module1_recurrence.features.scripts.thumb_tiles --root data/sj/tcga_other [--glob "*.svs"]
결과: {root}/tiles/{슬라이드}.csv
"""
import argparse
import math
from pathlib import Path

import numpy as np
import openslide
import pandas as pd
from PIL import Image, ImageFilter

from module1_recurrence.data.scripts.slide_qc import MAX_SIDE, measure
from module1_recurrence.features.scripts.tumor_tiles import PEN_DILATE, PEN_MAX, TILE_UM, to_grid

MIN_TISSUE = 0.5


def thumbnail(slide, qc_dir, name):
    path = qc_dir / "thumbs" / f"{name}.png"
    w, h = slide.dimensions
    down = 2 ** max(0, math.ceil(math.log2(max(w, h) / MAX_SIDE)))
    if path.exists():
        img = Image.open(path).convert("RGB")
        down = 2 ** round(math.log2(w / img.width))
    else:
        img = slide.get_thumbnail((math.ceil(w / down), math.ceil(h / down))).convert("RGB")
        path.parent.mkdir(parents=True, exist_ok=True)
        img.save(path)
    return img, down


def process(path, qc_dir, out):
    name = path.name.split(".")[0]
    slide = openslide.OpenSlide(str(path))
    mpp = float(slide.properties.get(openslide.PROPERTY_NAME_MPP_X))
    tile = round(TILE_UM / mpp)
    img, down = thumbnail(slide, qc_dir, name)
    _, pen = measure(img)
    hsv = np.asarray(img.convert("HSV")).astype(int)
    tissue_px = (hsv[..., 1] > 20) & (hsv[..., 2] < 235) & ~pen
    w, h = slide.dimensions
    nx, ny = math.ceil(w / tile), math.ceil(h / tile)
    pen_wide = np.asarray(Image.fromarray(pen.astype(np.uint8) * 255).filter(ImageFilter.MaxFilter(PEN_DILATE))) > 0
    frac, pen_frac = to_grid(tissue_px, nx, ny, down, tile), to_grid(pen_wide, nx, ny, down, tile)
    keep = (frac >= MIN_TISSUE) & (pen_frac <= PEN_MAX)
    ty, tx = np.nonzero(keep)
    df = pd.DataFrame({"x": tx * tile, "y": ty * tile, "tile": tile, "n_cells": 0, "tumor_ratio": np.nan,
                       "tissue_frac": frac[ty, tx].round(3), "pen_frac": pen_frac[ty, tx].round(3),
                       "tissue_tile": True, "tumor_tile": False})
    df.to_csv(out / f"{name}.csv", index=False)
    return name, len(df), mpp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True, help="slides/ · qc/ 가 있는 폴더")
    ap.add_argument("--glob", default="*.svs", help="슬라이드 파일 패턴 (HistologyHSI-BC는 '*.mrxs')")
    ap.add_argument("--slide_dir", default="slides")
    args = ap.parse_args()
    out = args.root / "tiles"
    out.mkdir(exist_ok=True)
    qc = pd.read_csv(args.root / "qc" / "slide_qc.csv").set_index("slide").stain if (args.root / "qc" / "slide_qc.csv").exists() else {}
    for p in sorted((args.root / args.slide_dir).glob(args.glob)):
        name = p.name.split(".")[0]
        if (out / f"{name}.csv").exists() or qc.get(name) == "IHC_suspect":
            continue
        try:
            print(*process(p, args.root / "qc", out), flush=True)
        except Exception as e:  # noqa: BLE001 (깨진 슬라이드는 건너뛰고 목록만 남김)
            print(name, "failed:", e, flush=True)


if __name__ == "__main__":
    main()
