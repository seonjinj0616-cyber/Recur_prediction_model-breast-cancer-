"""종양 타일 선택(tumor_tiles) 결과의 조직 타일을 ImageNet ResNet50으로 특징 벡터로 바꾼다 (B3).

인코더는 torchvision ImageNet ResNet50 (IMAGENET1K_V1)의 layer3까지 + 평균 풀링 → 타일당 1024차원.
타일은 tumor_tiles CSV와 같은 격자(128 µm, 원본 `tile` 픽셀)를 원본에서 읽어 224×224로 줄여 넣는다.
조직 타일 전부를 저장하고 `tumor_tile` 표시를 같이 남겨, 학습할 때 종양 타일 / 전체 조직 버전을 고른다.
슬라이드마다 파일 하나를 만들고, 환자 단위(DX1만 / 모든 슬라이드)로 묶는 것은 학습 단계에서 한다.
이미 .h5가 있는 슬라이드는 건너뛴다.
실행: python -m module1_recurrence.features.scripts.extract [--root DIR] [--device cuda] [--batch_size 256]
결과: {root}/features/{슬라이드}.h5
  features (N, 1024) float16, coords (N, 2) 타일 왼쪽 위 원본 픽셀 [x, y], tumor_tile (N,), tumor_ratio (N,),
  n_cells (N,), 속성 tile · mpp · encoder
"""
import argparse
import time
from pathlib import Path

import h5py
import numpy as np
import openslide
import pandas as pd
import torch
from PIL import Image
from torch.utils.data import DataLoader, Dataset
from torchvision import models, transforms

from module1_recurrence.data.scripts.gdc_fetch import RAW

ENCODER = "resnet50_imagenet1k_v1_layer3_avgpool"
PREP = transforms.Compose([
    transforms.ToTensor(),
    transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225]),
])


class Tiles(Dataset):
    """원본 슬라이드에서 타일을 읽는다. openslide 핸들은 작업자(worker)마다 따로 연다."""

    def __init__(self, svs, xy, tile):
        self.svs, self.xy, self.tile, self.slide = str(svs), xy, int(tile), None

    def __len__(self):
        return len(self.xy)

    def __getitem__(self, i):
        if self.slide is None:
            self.slide = openslide.OpenSlide(self.svs)
        x, y = self.xy[i]
        img = self.slide.read_region((int(x), int(y)), 0, (self.tile, self.tile)).convert("RGB")
        return PREP(img.resize((224, 224), Image.BILINEAR))


def encoder(device):
    m = models.resnet50(weights=models.ResNet50_Weights.IMAGENET1K_V1)
    net = torch.nn.Sequential(m.conv1, m.bn1, m.relu, m.maxpool, m.layer1, m.layer2, m.layer3,
                              torch.nn.AdaptiveAvgPool2d(1), torch.nn.Flatten())
    return net.eval().to(device)


@torch.no_grad()
def extract(net, svs, tiles, device, batch_size, workers):
    loader = DataLoader(Tiles(svs, tiles[["x", "y"]].to_numpy(), tiles.tile.iloc[0]),
                        batch_size=batch_size, num_workers=workers)
    out = []
    for batch in loader:
        with torch.autocast(device_type=device.split(":")[0], dtype=torch.float16, enabled=device != "cpu"):
            out.append(net(batch.to(device, non_blocking=True)).float().cpu().numpy())
    return np.concatenate(out).astype(np.float16)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=RAW.parent, help="slides/와 tiles/가 있는 폴더 (서버: data/sj/tcga)")
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--batch_size", type=int, default=256)
    ap.add_argument("--workers", type=int, default=8, help="타일 읽기 프로세스 수")
    args = ap.parse_args()

    out = args.root / "features"
    out.mkdir(parents=True, exist_ok=True)
    svs = {p.name.split(".")[0]: p for p in (args.root / "slides").glob("*.svs")}
    net = None

    for csv in sorted((args.root / "tiles").glob("TCGA-*.csv")):
        slide = csv.stem
        dest = out / f"{slide}.h5"
        if dest.exists() or slide not in svs:
            continue
        tiles = pd.read_csv(csv)
        tiles = tiles[tiles.tissue_tile].reset_index(drop=True)
        if tiles.empty:
            print(slide, "no tissue tiles")
            continue
        net = net or encoder(args.device)
        t0 = time.time()
        feats = extract(net, svs[slide], tiles, args.device, args.batch_size, args.workers)
        mpp = float(openslide.OpenSlide(str(svs[slide])).properties.get(openslide.PROPERTY_NAME_MPP_X, "nan"))

        tmp = dest.with_suffix(".part")
        with h5py.File(tmp, "w") as f:
            f["features"] = feats
            f["coords"] = tiles[["x", "y"]].to_numpy(dtype=np.int64)
            f["tumor_tile"] = tiles.tumor_tile.to_numpy()
            f["tumor_ratio"] = tiles.tumor_ratio.to_numpy(dtype=np.float32)
            f["n_cells"] = tiles.n_cells.to_numpy(dtype=np.int32)
            f.attrs.update(tile=int(tiles.tile.iloc[0]), mpp=mpp, encoder=ENCODER)
        tmp.rename(dest)  # 다 쓴 파일만 .h5 이름을 가진다
        print(f"{slide} tiles {len(tiles)} (tumor {int(tiles.tumor_tile.sum())}) {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
