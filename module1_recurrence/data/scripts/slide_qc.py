"""TCGA-BRCA 진단 슬라이드 전체를 원본 없이 썸네일로 점검한다 (염색 종류 · 펜 표시).

GDC 슬라이드 뷰어용 타일 API(DeepZoom)에서 긴 변 1024픽셀 이하 썸네일만 받는다 (슬라이드당 약 1 MB).
- 염색: 조직 픽셀 중 H&E 색(보라~분홍) 비율이 낮고 갈색(DAB)이 있으면 IHC 의심 → 눈으로 확인한 뒤 분석에서 뺀다.
  H&E 색 비율만 낮으면 에오신이 바랜 H&E (`H&E_faded`)로 따로 표시한다.
- 펜: 초록 · 파랑 · 검정 펜 픽셀 비율과 마스크(png)를 남긴다 → 타일 단계에서 펜이 지나가는 타일을 뺀다.
실행: python -m module1_recurrence.data.scripts.slide_qc [--limit N] [--root DIR]
      GDC 썸네일이 안 된 슬라이드: ... slide_qc --root DIR --from_svs TCGA-AC-A6NO-01Z-00-DX1 (받아 둔 원본으로 다시 판정)
결과: {root}/qc/slide_qc.csv, {root}/qc/thumbs/*.png, {root}/qc/pen_masks/*.png,
      {root}/qc/{IHC_suspect,HE_faded,check_thumb,pen,red_outline,red_remnant,red_tissue_red}.png (종류별 썸네일 모음)
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
from scipy.ndimage import label

from module1_recurrence.data.scripts.gdc_fetch import API, RAW

MAX_SIDE = 1024
RETRIES = 3
LOW_HE = 0.5  # 조직 중 H&E 색(보라~분홍) 비율이 이보다 낮으면 청회색 슬라이드
IHC_BROWN = 0.02  # 청회색이면서 갈색(DAB)이 이만큼 있으면 IHC 의심. 1,133장 중 갈색 있는 것은 5L-AAT0(0.078)뿐,
# 나머지 청회색 슬라이드(대부분 A8 기관)는 에오신이 바랜 H&E였음
MIN_TISSUE = 0.02  # 썸네일 조직 비율이 이보다 낮으면 썸네일이 깨졌거나 조직이 거의 없음 → 원본으로 다시 확인


def _get(url):
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read()
        except OSError as e:
            # 400 같은 요청 오류는 다시 해도 같으므로 바로 넘긴다
            if attempt == RETRIES - 1 or (isinstance(e, urllib.error.HTTPError) and e.code < 500):
                raise
            time.sleep(5 * (attempt + 1))


def pyramid(file_id):
    """DeepZoom 메타데이터 → 원본 크기, 타일 크기 · 겹침, 그리고 긴 변이 MAX_SIDE 이하가 되는 축소 배율."""
    meta = json.loads(_get(f"{API}/tile/metadata/{file_id}"))
    w, h, ts, ov = int(meta["Width"]), int(meta["Height"]), int(meta["TileSize"]), int(meta["Overlap"])
    down = 2 ** max(0, math.ceil(math.log2(max(w, h) / MAX_SIDE)))
    return w, h, ts, ov, down


def thumbnail(file_id):
    """긴 변이 MAX_SIDE 이하인 단계부터 받아 보고, GDC가 그 단계를 못 그리면(400) 한 단계씩 낮춘다."""
    w, h, ts, ov, down = pyramid(file_id)
    for _ in range(3):
        try:
            return _stitch(file_id, w, h, ts, ov, down), down
        except urllib.error.HTTPError as e:
            if e.code != 400:
                raise
            down *= 2
    raise OSError(f"no renderable level down to 1/{down}")


def _stitch(file_id, w, h, ts, ov, down):
    level = math.ceil(math.log2(max(w, h))) - int(math.log2(down))  # 첫 항은 원본 해상도 단계
    lw, lh = math.ceil(w / down), math.ceil(h / down)
    img = Image.new("RGB", (lw, lh), "white")
    for ty in range(math.ceil(lh / ts)):
        for tx in range(math.ceil(lw / ts)):
            tile = Image.open(io.BytesIO(_get(f"{API}/tile/{file_id}?level={level}&x={tx}&y={ty}"))).convert("RGB")
            # 첫 행 · 열이 아니면 타일 앞쪽에 겹침(overlap) 픽셀이 붙어 있다
            img.paste(tile, (tx * ts - (ov if tx else 0), ty * ts - (ov if ty else 0)))
    return img


def measure(img):
    """HSV 색 범위로 조직 · 분홍 · 갈색 · 펜 픽셀을 센다. 비율은 조직(펜 제외) 픽셀 기준."""
    hsv = np.asarray(img.convert("HSV")).astype(int)
    hue, sat, val = hsv[..., 0] * 360 // 255, hsv[..., 1], hsv[..., 2]
    green = (hue >= 90) & (hue <= 200) & (sat > 60)
    blue = (hue >= 200) & (hue <= 260) & (sat > 110) & (val < 200)
    dark = val < 60
    # 진한 빨강: 빨간 펜 잉크와 혈액 · 진한 에오신이 함께 걸린다 → red_kind로 나눈다
    red = ((hue >= 345) | (hue <= 10)) & (sat > 170) & (val > 80)
    kind = red_kind(red, sat)
    pen = green | blue | dark | (red if kind in ("outline", "remnant") else False)
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
        "red_pen_frac": red.mean(),
        "red_kind": kind,
        "dark_frac": dark.mean(),
    }
    return stats, pen


def red_kind(red, sat):
    """진한 빨강 픽셀이 무엇인지 나눈다. 빨강이 0.5% 이상인 39장을 눈으로 나눈 결과와 39장 모두 일치하도록 정한 기준.

    - tissue_red: 조직 자체가 빨강 (혈액 · 진한 에오신). 펜 잉크만큼 선명한 픽셀(채도 > 220)이 거의 없음 (≤ 0.02, 펜은 ≥ 0.23)
    - broken: 빨간 덩어리 하나가 외곽 상자를 꽉 채움 → 썸네일 저장 이상 (TCGA-OL-A5RZ)
    - outline: 빨강 5% 이상 → 종양을 둘러싼 빨간 펜 선 (6.4–16.7%, 대부분 E9 기관)
    - remnant: 그 밖 → 지워진 펜 흔적 (≤ 3.9%). 잉크 픽셀만 빼고 영역 표시로는 쓰지 않음
    """
    frac = red.mean()
    if frac < 0.005:
        return "none"
    if (red & (sat > 220)).sum() / red.sum() < 0.1:
        return "tissue_red"
    lab, n = label(red)
    sizes = np.bincount(lab.ravel())[1:]
    ys, xs = np.nonzero(lab == sizes.argmax() + 1)
    if sizes.max() / red.sum() > 0.9 and sizes.max() / ((np.ptp(xs) + 1) * (np.ptp(ys) + 1)) > 0.9:
        return "broken"
    return "outline" if frac >= 0.05 else "remnant"


def run_one(row, out):
    name = row.file_name.split(".")[0]  # 예: TCGA-A1-A0SB-01Z-00-DX1
    thumb_path = out / "thumbs" / f"{name}.png"
    try:
        if thumb_path.exists():  # 다시 돌릴 때는 썸네일 대신 메타데이터만 받아 축소 배율을 되찾는다
            img = Image.open(thumb_path).convert("RGB")
            down = 2 ** round(math.log2(pyramid(row.file_id)[0] / img.width))
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


def from_svs(names, root):
    """GDC 썸네일이 깨졌거나(check_thumb) 받지 못한(error) 슬라이드를 받아 둔 원본 svs로 다시 판정해 qc 표를 고친다."""
    import openslide
    out = root / "qc"
    qc = pd.read_csv(out / "slide_qc.csv")
    for name in names:
        o = openslide.OpenSlide(str(next((root / "slides").glob(f"{name}.*.svs"))))
        w, h = o.dimensions
        down = 2 ** max(0, math.ceil(math.log2(max(w, h) / MAX_SIDE)))
        img = o.get_thumbnail((math.ceil(w / down), math.ceil(h / down))).convert("RGB")
        img.save(out / "thumbs" / f"{name}.png")
        stats, pen = measure(img)
        Image.fromarray(pen.astype(np.uint8) * 255).save(out / "pen_masks" / f"{name}.png")
        i = qc.index[qc.slide == name][0]
        for k, v in {"thumb_w": img.width, "thumb_h": img.height, "downsample": down, **stats}.items():
            qc.loc[i, k] = v
        qc.loc[i, "error"] = None
        low_he = qc.loc[i, "he_frac"] < LOW_HE
        stain = "IHC_suspect" if low_he and qc.loc[i, "brown_frac"] >= IHC_BROWN else "H&E_faded" if low_he else "H&E"
        qc.loc[i, "stain"] = "check_thumb" if stats["tissue_frac"] < MIN_TISSUE or stats["red_kind"] == "broken" else stain
        red_ink = stats["red_kind"] in ("outline", "remnant")
        qc.loc[i, "pen"] = (stats["green_pen_frac"] + stats["blue_pen_frac"] + (stats["red_pen_frac"] if red_ink else 0)) > 0.005
        print(qc.loc[i, ["slide", "stain", "tissue_frac", "he_frac", "brown_frac", "pen", "downsample"]].to_dict())
    qc.to_csv(out / "slide_qc.csv", index=False)


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
    ap.add_argument("--from_svs", nargs="+", metavar="SLIDE", help="이 슬라이드만 받아 둔 원본 svs로 다시 판정 (예: TCGA-AC-A6NO-01Z-00-DX1)")
    args = ap.parse_args()
    if args.from_svs:
        return from_svs(args.from_svs, args.root)

    out = args.root / "qc"
    for d in ("thumbs", "pen_masks"):
        (out / d).mkdir(parents=True, exist_ok=True)

    # 환자 단위로 거르지 않고 TCGA-BRCA 진단 슬라이드 전부 (1,047명 코호트를 포함)
    dx = pd.read_csv(args.root / "raw" / "dx_slides.tsv", sep="\t").rename(
        columns={"cases.0.submitter_id": "patient_id"}).sort_values("file_name").head(args.limit)

    with ThreadPoolExecutor(max_workers=8) as ex:
        rows = list(ex.map(lambda r: run_one(r, out), dx.itertuples()))

    qc = pd.DataFrame(rows)
    if "error" not in qc:
        qc["error"] = None
    low_he = qc.he_frac < LOW_HE
    qc["stain"] = "H&E"
    qc.loc[low_he, "stain"] = "H&E_faded"
    qc.loc[low_he & (qc.brown_frac >= IHC_BROWN), "stain"] = "IHC_suspect"
    qc.loc[(qc.tissue_frac < MIN_TISSUE) | (qc.red_kind == "broken"), "stain"] = "check_thumb"
    qc.loc[qc.error.notna(), "stain"] = "error"
    red_ink = qc.red_kind.isin(["outline", "remnant"])
    qc["pen"] = (qc.green_pen_frac + qc.blue_pen_frac + qc.red_pen_frac.where(red_ink, 0)) > 0.005
    qc.to_csv(out / "slide_qc.csv", index=False)

    for kind in ("IHC_suspect", "H&E_faded", "check_thumb"):
        s = qc[qc.stain == kind].sort_values("he_frac")
        contact_sheet([out / "thumbs" / f"{n}.png" for n in s.slide], out / f"{kind.replace('&', '')}.png")
    contact_sheet([out / "thumbs" / f"{n}.png" for n in qc[qc.pen].slide], out / "pen.png")
    for kind in ("outline", "remnant", "tissue_red"):
        s = qc[qc.red_kind == kind].sort_values("red_pen_frac", ascending=False)
        contact_sheet([out / "thumbs" / f"{n}.png" for n in s.slide], out / f"red_{kind}.png")
    print(qc.stain.value_counts().to_dict(), "| pen (> 0.5%)", int(qc.pen.sum()),
          "| red", qc.red_kind.value_counts().drop("none", errors="ignore").to_dict())
    print(qc[qc.stain.isin(["IHC_suspect", "check_thumb", "error"])][["slide", "stain", "he_frac", "brown_frac"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
