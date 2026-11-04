"""B1–B3 묶음 실행: 슬라이드 받기 → CellViT 세포 분류 → 종양 타일 선택 → 특징 추출.

슬라이드는 지우지 않는다 (종양 범위 지정 방법 비교가 끝난 뒤 정함). 각 단계는 이미 끝난 슬라이드를 건너뛰므로
같은 명령을 다시 실행하면 이어서 진행한다. team_project 환경에서 실행하고, CellViT는 cellvit 환경의 실행 파일을 부른다.
썸네일 QC가 없거나 확인 필요(check_thumb · error)인 슬라이드는 CellViT를 돌리지 않고 목록만 출력한다.
CellViT 결과 중 핵 윤곽 파일(cells.json · cells.geojson, 슬라이드당 수백 MB)은 지우고
핵 중심 파일(cell_detection.json · .geojson, QuPath 확인용)만 남긴다 (--keep_contours로 유지).
실행: python -m module1_recurrence.features.scripts.run_batch --root data/sj/tcga --batch 20 --gpu 1 [--all_slides]
결과: {root}/slides, {root}/cellvit, {root}/tiles, {root}/features, {root}/logs/{슬라이드}.log
"""
import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

from module1_recurrence.data.scripts.gdc_fetch import RAW

CELLVIT = Path.home() / "miniconda3" / "envs" / "cellvit" / "bin" / "cellvit-inference"


def step(module, *args, env=None):
    """같은 파이썬으로 다른 스크립트를 실행한다. 실패하면 멈춘다."""
    subprocess.run([sys.executable, "-m", f"module1_recurrence.{module}", *map(str, args)], check=True, env=env)


def run_cellvit(svs, out_dir, log_dir, cellvit, env, keep_contours):
    t0 = time.time()
    with open(log_dir / f"{svs.name.split('.')[0]}.log", "w") as log:
        r = subprocess.run([str(cellvit), "--model", "SAM", "--nuclei_taxonomy", "pannuke", "--outdir", str(out_dir),
                            "--geojson", "process_wsi", "--wsi_path", str(svs)], stdout=log, stderr=subprocess.STDOUT, env=env)
    res = out_dir / svs.stem
    if r.returncode == 0 and not keep_contours:
        for f in ("cells.json", "cells.geojson"):
            (res / f).unlink(missing_ok=True)
    return r.returncode, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=RAW.parent, help="데이터 폴더 (서버: data/sj/tcga)")
    ap.add_argument("--batch", type=int, default=20, help="이번에 새로 받을 슬라이드 수 (0이면 받지 않고 남은 처리만)")
    ap.add_argument("--gpu", default="0", help="CellViT · 특징 추출에 쓸 GPU 번호")
    ap.add_argument("--all_slides", action="store_true", help="환자의 DX 슬라이드 전부 (기본: 첫 번째만)")
    ap.add_argument("--cellvit", type=Path, default=CELLVIT, help="cellvit-inference 실행 파일")
    ap.add_argument("--keep_contours", action="store_true", help="CellViT 핵 윤곽 파일을 지우지 않음")
    args = ap.parse_args()

    root = args.root
    env = {**os.environ, "CUDA_VISIBLE_DEVICES": args.gpu}
    cell_dir, log_dir = root / "cellvit", root / "logs"
    cell_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    if args.batch:
        step("data.scripts.slide_download", "--root", root, "--batch", args.batch, *(["--all_slides"] if args.all_slides else []))

    qc = pd.read_csv(root / "qc" / "slide_qc.csv").set_index("slide").stain
    done = {d.name.split(".")[0] for d in cell_dir.iterdir() if (d / "cell_detection.json").exists()}
    held = []
    for svs in sorted((root / "slides").glob("*.svs")):
        slide = svs.name.split(".")[0]
        if slide in done:
            continue
        stain = qc.get(slide, "no_qc")
        if stain in ("IHC_suspect", "check_thumb", "error", "no_qc"):
            held.append((slide, stain))
            continue
        code, sec = run_cellvit(svs, cell_dir, log_dir, args.cellvit, env, args.keep_contours)
        print(f"cellvit {slide} {'ok' if code == 0 else f'FAILED ({code}), see logs/{slide}.log'} {sec / 60:.1f} min",
              flush=True)

    step("features.scripts.tumor_tiles", "--root", root)
    step("features.scripts.extract", "--root", root, env=env)
    if held:
        print("held (원본 확인 필요):", held)


if __name__ == "__main__":
    main()
