"""CellViT 세포 분류 결과로 슬라이드를 타일로 나눠 종양 핵 비율을 세고, 종양 타일을 고른다 (B2, 단계 1).

타일은 20배 256×256에 해당하는 128 µm 정사각형. 원본 픽셀로는 슬라이드 mpp에 따라 다르다 (40배 0.25 mpp → 약 508픽셀).
원점 (0, 0)부터 겹침 없이 나누고, 크기를 CSV `tile` 열에 남긴다 (B3 extract.py가 같은 격자를 쓴다).
펜이 지나가는 타일은 빼고(옅은 펜도 분류를 왜곡해서 마스크를 넓혀 씀), IHC 의심 슬라이드는 건너뛴다.
먼저 data.scripts.slide_qc로 {root}/qc/를 만들고, CellViT 결과를 {root}/cellvit/{슬라이드}/에 둔다.
이미 CSV가 있는 슬라이드는 건너뛴다 (--redo로 다시 계산).
실행: python -m module1_recurrence.features.scripts.tumor_tiles [--root DIR] [--min_ratio 0.5] [--min_cells 20] [--redo]
결과: {root}/tiles/{슬라이드}.csv (타일별 핵 수), {root}/tiles/summary.csv,
      {root}/tiles/{슬라이드}.jpg (왼쪽 썸네일 · 오른쪽 선택 결과: 빨강 종양, 파랑 그 밖 조직, 회색 펜 제외)
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFilter
from scipy.ndimage import label
from scipy.spatial import ConvexHull

from module1_recurrence.data.scripts.gdc_fetch import RAW

TILE_UM = 128  # 20배(0.5 mpp) 256픽셀
TYPES = {1: "neoplastic", 2: "inflammatory", 3: "connective", 4: "dead", 5: "epithelial"}
PEN_DILATE = 5  # 썸네일 픽셀 (홀수). 썸네일 1픽셀 = 원본 64–128픽셀
PEN_MAX = 0.01  # 타일 중 펜 비율이 이보다 크면 제외
OUTLINE_MIN = 0.05  # 펜 덩어리의 외곽 상자가 썸네일의 이 비율보다 커야 "원"으로 봄
COLORS = {"tumor": (220, 30, 30), "tissue": (40, 90, 230), "pen": (130, 130, 130)}


def load_cells(path):
    """cell_detection.json → (x, y, type) 배열과 원본 mpp. centroid는 [x, y] (bbox는 [y, x] 순서라 쓰지 않음)."""
    det = json.loads(path.read_text())
    xy = np.array([c["centroid"] for c in det["cells"]], dtype=np.int64)
    t = np.array([c["type"] for c in det["cells"]], dtype=np.int64)
    return xy[:, 0], xy[:, 1], t, det["wsi_metadata"]["base_mpp"]


def count_tiles(x, y, t, nx, ny, tile):
    """타일 (ty, tx)마다 종류별 핵 수. counts[..., k]는 type k (0은 비워 둠)."""
    tx, ty = np.clip(x // tile, 0, nx - 1), np.clip(y // tile, 0, ny - 1)
    counts = np.zeros((ny, nx, 6), dtype=np.int32)
    np.add.at(counts, (ty, tx, t), 1)
    return counts


def to_grid(mask, nx, ny, down, tile):
    """썸네일 크기 마스크를 타일 격자로 줄여 타일별 비율을 낸다."""
    # 썸네일은 원본 / down 크기라 격자 끝과 맞추려면 먼저 격자 크기 × (tile / down)로 맞춘다
    scale = tile / down
    w, h = max(1, round(nx * scale)), max(1, round(ny * scale))
    m = Image.fromarray(mask.astype(np.uint8) * 255).crop((0, 0, w, h))
    return np.asarray(m.resize((nx, ny), Image.BOX), dtype=float) / 255


def outline_inside(pen):
    """초록 원 같은 펜 선의 안쪽 (보조 근거로만 씀).

    손으로 그린 원은 끝이 벌어져 있어 구멍 채우기가 안 되므로, 큰 펜 덩어리마다 볼록 껍질(convex hull)을 칠한다.
    짧은 선 · 점 표시는 껍질이 작아서 빠진다.
    """
    lab, n = label(pen)
    inside = Image.new("L", (pen.shape[1], pen.shape[0]), 0)
    draw = ImageDraw.Draw(inside)
    for k in range(1, n + 1):
        ys, xs = np.nonzero(lab == k)
        if len(xs) < 3 or (np.ptp(xs) + 1) * (np.ptp(ys) + 1) < OUTLINE_MIN * pen.size:
            continue
        pts = np.column_stack([xs, ys])
        hull = pts[ConvexHull(pts).vertices]
        draw.polygon([tuple(p) for p in hull], fill=255)
    return (np.asarray(inside) > 0) & ~pen


def overlay(thumb, labels, nx, ny, down, tile, dest):
    """썸네일 위에 타일 판정을 반투명으로 칠한다."""
    rgb = np.zeros((ny, nx, 3), dtype=np.uint8)
    alpha = np.zeros((ny, nx), dtype=np.uint8)
    for name, color in COLORS.items():
        rgb[labels == name], alpha[labels == name] = color, 110
    scale = tile / down
    size = (round(nx * scale), round(ny * scale))
    layer = Image.fromarray(rgb).resize(size, Image.NEAREST)
    a = Image.fromarray(alpha).resize(size, Image.NEAREST)
    base = thumb.copy()
    base.paste(layer, (0, 0), a)
    # 왼쪽 원본 썸네일, 오른쪽 선택 결과 (눈으로 비교용)
    pair = Image.new("RGB", (thumb.width * 2 + 10, thumb.height), "white")
    pair.paste(thumb, (0, 0))
    pair.paste(base, (thumb.width + 10, 0))
    pair.save(dest, quality=90)


def process(slide, det_path, qc_row, qc_dir, out, min_ratio, min_cells):
    t0 = time.time()
    down = qc_row.downsample
    thumb = Image.open(qc_dir / "thumbs" / f"{slide}.png").convert("RGB")
    x, y, t, mpp = load_cells(det_path)
    tile = round(TILE_UM / mpp)
    nx, ny = int(np.ceil(thumb.width * down / tile)), int(np.ceil(thumb.height * down / tile))
    counts = count_tiles(x, y, t, nx, ny, tile)

    pen = np.asarray(Image.open(qc_dir / "pen_masks" / f"{slide}.png")) > 0
    pen_wide = np.asarray(Image.fromarray(pen.astype(np.uint8) * 255).filter(ImageFilter.MaxFilter(PEN_DILATE))) > 0
    pen_frac = to_grid(pen_wide, nx, ny, down, tile)
    inside = to_grid(outline_inside(pen_wide), nx, ny, down, tile) > 0.5

    n = counts[..., 1:].sum(-1)
    ratio = np.divide(counts[..., 1], n, out=np.zeros(n.shape), where=n > 0)
    is_pen = pen_frac > PEN_MAX
    tissue = (n >= min_cells) & ~is_pen
    tumor = tissue & (ratio >= min_ratio)

    labels = np.full((ny, nx), "", dtype=object)
    labels[tissue], labels[tumor], labels[is_pen & (n > 0)] = "tissue", "tumor", "pen"
    overlay(thumb, labels, nx, ny, down, tile, out / f"{slide}.jpg")

    ty, tx = np.nonzero(n > 0)
    df = pd.DataFrame({"x": tx * tile, "y": ty * tile, "tile": tile, "n_cells": n[ty, tx],
                       **{f"n_{v}": counts[ty, tx, k] for k, v in TYPES.items()},
                       "tumor_ratio": ratio[ty, tx].round(4), "pen_frac": pen_frac[ty, tx].round(4),
                       "in_pen_outline": inside[ty, tx], "tissue_tile": tissue[ty, tx], "tumor_tile": tumor[ty, tx]})
    df.to_csv(out / f"{slide}.csv", index=False)

    has_outline = bool(inside.any())
    return {"slide": slide, "stain": qc_row.stain, "mpp": mpp, "tile": tile, "cells": len(t), "tiles_with_cells": len(df),
            "tissue_tiles": int(tissue.sum()), "tumor_tiles": int(tumor.sum()), "pen_tiles": int((is_pen & (n > 0)).sum()),
            "tumor_in_outline": round(float((tumor & inside).sum() / max(tumor.sum(), 1)), 3) if has_outline else None,
            "seconds": round(time.time() - t0, 1)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=RAW.parent, help="qc/와 cellvit_test/가 있는 폴더 (서버: data/sj/tcga)")
    ap.add_argument("--min_ratio", type=float, default=0.5, help="종양 타일 기준: Neoplastic / 전체 핵")
    ap.add_argument("--min_cells", type=int, default=20, help="조직 타일 기준: 타일 안 최소 핵 수")
    ap.add_argument("--redo", action="store_true", help="CSV가 있어도 다시 계산 (기준값을 바꿨을 때)")
    args = ap.parse_args()

    qc_dir, out = args.root / "qc", args.root / "tiles"
    out.mkdir(parents=True, exist_ok=True)
    qc = pd.read_csv(qc_dir / "slide_qc.csv").set_index("slide")
    summary_path = out / "summary.csv"
    old = pd.read_csv(summary_path) if summary_path.exists() and not args.redo else pd.DataFrame(columns=["slide"])

    rows = []
    for det in sorted((args.root / "cellvit").glob("*/cell_detection.json")):
        slide = det.parent.name.split(".")[0]  # 예: TCGA-A1-A0SB-01Z-00-DX1
        if not args.redo and slide in set(old.slide):
            continue
        row = qc.loc[slide]
        if row.stain == "IHC_suspect":
            rows.append({"slide": slide, "stain": row.stain})
            print(slide, "skip (IHC suspect)")
            continue
        r = process(slide, det, row, qc_dir, out, args.min_ratio, args.min_cells)
        rows.append(r)
        print(r)

    pd.concat([old, pd.DataFrame(rows)]).to_csv(summary_path, index=False)


if __name__ == "__main__":
    main()
