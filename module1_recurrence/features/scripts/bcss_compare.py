"""같은 타일 격자에서 CellViT 종양 타일과 BCSS 조직 분할(TIAToolbox fcn_resnet50_unet-bcss)의 종양 판정을 비교한다.

BCSS 모델은 0.25 mpp 1024×1024를 넣으면 가운데 512×512 (128 µm)의 조직 종류를 낸다 → tumor_tiles 타일(128 µm)과
크기가 같다. 타일마다 중심을 맞춰 256 µm를 읽어 1024로 줄이고, 가운데 출력에서 Tumour 픽셀 비율을 센다.
tiatoolbox 환경에서 실행한다 (먼저 tumor_tiles가 {root}/tiles/{슬라이드}.csv를 만들어야 함). 끝난 슬라이드는 건너뛴다.
실행: python -m module1_recurrence.features.scripts.bcss_compare --root data/sj/tcga [--slides A B ...] [--min_frac 0.5]
결과: {root}/tiles/bcss/{슬라이드}.csv (타일별 BCSS 5종 비율), {root}/tiles/bcss/{슬라이드}.jpg
      (썸네일 | CellViT 종양 | BCSS 종양 | 일치: 빨강 둘 다, 주황 CellViT만, 보라 BCSS만), {root}/tiles/bcss/summary.csv
"""
import argparse
import time
from pathlib import Path

import numpy as np
import openslide
import pandas as pd
import torch
from PIL import Image
from torch.utils.data import DataLoader, Dataset

CLASSES = ["tumour", "stroma", "inflammatory", "necrosis", "others"]  # pretrained_model.yaml class_dict 순서
MODEL_MPP, IN_PX = 0.25, 1024


class Windows(Dataset):
    """타일 중심에 맞춘 256 µm 창을 읽어 1024×1024 (0.25 mpp)로 맞춘다. 모델 입력은 NHWC uint8."""

    def __init__(self, svs, tiles, mpp):
        self.svs, self.slide = str(svs), None
        self.size = round(IN_PX * MODEL_MPP / mpp)
        self.xy = tiles[["x", "y"]].to_numpy() + tiles.tile.to_numpy()[:, None] // 2 - self.size // 2

    def __len__(self):
        return len(self.xy)

    def __getitem__(self, i):
        if self.slide is None:
            self.slide = openslide.OpenSlide(self.svs)
        x, y = self.xy[i]
        img = self.slide.read_region((int(x), int(y)), 0, (self.size, self.size)).convert("RGB")
        return torch.from_numpy(np.asarray(img.resize((IN_PX, IN_PX), Image.BILINEAR)).copy())


def panel(thumb, mask, color, down, tile):
    """썸네일 위에 타일 마스크를 칠한 그림 하나."""
    ny, nx = mask.shape[:2]
    size = (round(nx * tile / down), round(ny * tile / down))
    rgb = Image.fromarray(np.where(mask[..., None], np.array(color, np.uint8), 0).astype(np.uint8)).resize(size, Image.NEAREST)
    a = Image.fromarray((mask.any(-1) if mask.ndim == 3 else mask).astype(np.uint8) * 120).resize(size, Image.NEAREST)
    base = thumb.copy()
    base.paste(rgb, (0, 0), a)
    return base


def figure(thumb, df, down, dest):
    tile = int(df.tile.iloc[0])
    nx, ny = int(np.ceil(thumb.width * down / tile)), int(np.ceil(thumb.height * down / tile))
    gx, gy = (df.x // tile).to_numpy(), (df.y // tile).to_numpy()
    grids = {}
    for name in ("cv", "bc"):
        g = np.zeros((ny, nx), bool)
        g[gy, gx] = df[f"{name}_tumor"].to_numpy()
        grids[name] = g
    both, cv_only, bc_only = grids["cv"] & grids["bc"], grids["cv"] & ~grids["bc"], ~grids["cv"] & grids["bc"]
    agree = np.zeros((ny, nx, 3), np.uint8)
    agree[both], agree[cv_only], agree[bc_only] = (220, 30, 30), (255, 150, 0), (140, 40, 200)
    imgs = [thumb, panel(thumb, grids["cv"], (220, 30, 30), down, tile), panel(thumb, grids["bc"], (220, 30, 30), down, tile)]
    a_img = thumb.copy()
    size = (round(nx * tile / down), round(ny * tile / down))
    a_img.paste(Image.fromarray(agree).resize(size, Image.NEAREST), (0, 0),
                Image.fromarray(((both | cv_only | bc_only) * 140).astype(np.uint8)).resize(size, Image.NEAREST))
    imgs.append(a_img)
    out = Image.new("RGB", (thumb.width * 4 + 30, thumb.height), "white")
    for i, im in enumerate(imgs):
        out.paste(im, (i * (thumb.width + 10), 0))
    out.save(dest, quality=90)


def kappa(a, b):
    po = (a == b).mean()
    pe = a.mean() * b.mean() + (1 - a.mean()) * (1 - b.mean())
    return (po - pe) / (1 - pe) if pe < 1 else 1.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True, help="slides/, tiles/, qc/가 있는 폴더 (서버: data/sj/tcga)")
    ap.add_argument("--slides", nargs="*", help="비교할 슬라이드 (예: TCGA-A1-A0SB-01Z-00-DX1). 없으면 tiles/ 전부")
    ap.add_argument("--min_frac", type=float, default=0.5, help="BCSS 종양 타일 기준: 타일 중 Tumour 픽셀 비율")
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--batch_size", type=int, default=8)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    out = args.root / "tiles" / "bcss"
    out.mkdir(parents=True, exist_ok=True)
    svs = {p.name.split(".")[0]: p for p in (args.root / "slides").glob("*.svs")}
    qc = pd.read_csv(args.root / "qc" / "slide_qc.csv").set_index("slide")
    names = args.slides or [p.stem for p in sorted((args.root / "tiles").glob("TCGA-*.csv"))]
    from tiatoolbox.models.architecture import get_pretrained_model  # tiatoolbox 환경에만 있음
    model, _ = get_pretrained_model("fcn_resnet50_unet-bcss")
    model = model.to(args.device).eval()

    rows = []
    for slide in names:
        dest = out / f"{slide}.csv"
        if slide not in svs or not (args.root / "tiles" / f"{slide}.csv").exists():
            print(slide, "skip (no slide or tiles)")
            continue
        df = pd.read_csv(args.root / "tiles" / f"{slide}.csv")
        df = df[df.tissue_tile].reset_index(drop=True)
        if not dest.exists():
            t0 = time.time()
            mpp = float(openslide.OpenSlide(str(svs[slide])).properties[openslide.PROPERTY_NAME_MPP_X])
            loader = DataLoader(Windows(svs[slide], df, mpp), batch_size=args.batch_size, num_workers=args.workers)
            fracs = []
            for batch in loader:
                cls = model.infer_batch(model, batch, device=args.device).argmax(-1)  # (N, 512, 512)
                fracs.append(np.stack([(cls == k).mean((1, 2)) for k in range(len(CLASSES))], 1))
            df[[f"bcss_{c}" for c in CLASSES]] = np.concatenate(fracs).round(4)
            df.to_csv(dest, index=False)
            print(f"{slide} bcss {len(df)} tiles {time.time() - t0:.0f}s", flush=True)
        df = pd.read_csv(dest)
        df["cv_tumor"], df["bc_tumor"] = df.tumor_tile, df.bcss_tumour >= args.min_frac
        figure(Image.open(args.root / "qc" / "thumbs" / f"{slide}.png").convert("RGB"), df,
               qc.loc[slide, "downsample"], out / f"{slide}.jpg")
        cv, bc = df.cv_tumor.to_numpy(), df.bc_tumor.to_numpy()
        rows.append({"slide": slide, "stain": qc.loc[slide, "stain"], "tissue_tiles": len(df),
                     "cellvit_tumor": int(cv.sum()), "bcss_tumor": int(bc.sum()),
                     "both": int((cv & bc).sum()), "agree": round(float((cv == bc).mean()), 3),
                     "kappa": round(float(kappa(cv, bc)), 3),
                     "dice": round(float(2 * (cv & bc).sum() / max(cv.sum() + bc.sum(), 1)), 3),
                     "spearman": round(float(df.tumor_ratio.corr(df.bcss_tumour, method="spearman")), 3)})
        print(rows[-1], flush=True)

    pd.DataFrame(rows).to_csv(out / "summary.csv", index=False)


if __name__ == "__main__":
    main()
