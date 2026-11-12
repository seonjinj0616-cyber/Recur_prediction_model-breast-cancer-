"""다른 작업이 이미 받은 같은 GDC 슬라이드를 다시 받지 않고 바로가기(심볼릭 링크)로 쓴다.

TCGA 다른 아형 532명(검증용)과 모듈 4-B(HER2 양성 · 삼중음성)가 같은 DX1 슬라이드를 받는다 (2026-11-11 확인: 235장 겹침).
겹치는 file_id 목록을 {root}/exclude_shared.txt로 쓰고(slide_download --exclude에 넘김),
상대 폴더에 크기까지 다 받은 파일이 있으면 {root}/slides/에 링크를 만든다. 상대 쪽 파일은 읽기만 한다.
상대 쪽이 나중에 슬라이드를 지우면 링크가 끊기므로, 특징 추출이 끝날 때까지는 지우지 않아야 한다.
실행: python -m module1_recurrence.data.scripts.link_shared_slides --root data/sj/tcga_other
      --manifest ~/module4b/cohort/gdc_manifest_her2_tnbc.txt --shared_dir ~/tcga4b/hopt/slides
"""
import argparse
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True, help="상대 작업의 GDC manifest (id · filename · size)")
    ap.add_argument("--shared_dir", type=Path, required=True, help="상대 작업이 슬라이드를 받는 폴더")
    ap.add_argument("--replace_own", action="store_true",
                    help="우리가 이미 받은 겹치는 슬라이드도 지우고 상대 파일 링크로 바꿈 (상대 파일이 크기까지 다 받아진 경우만)")
    args = ap.parse_args()

    ours = pd.read_csv(args.root / "manifest.csv")
    ours_ids = set(ours.slide_file_ids.str.split(";").str[0])
    other = pd.read_csv(args.manifest, sep="\t")
    shared = other[other.id.isin(ours_ids)]
    (args.root / "exclude_shared.txt").write_text("\n".join(shared.id) + "\n")

    slide_dir = args.root / "slides"
    linked = waiting = mine = freed = 0
    for name, size in zip(shared.filename, shared["size"]):
        dest, src = slide_dir / name, args.shared_dir / name
        src_ok = src.exists() and src.stat().st_size == size
        if dest.exists() and not dest.is_symlink() and args.replace_own and src_ok:
            freed += dest.stat().st_size
            dest.unlink()
            dest.symlink_to(src)
            linked += 1
        elif dest.exists() and not dest.is_symlink():
            mine += 1  # 이미 우리가 받은 것은 그대로 둔다
        elif src_ok:
            if not dest.is_symlink():
                dest.symlink_to(src)
            linked += 1
        else:
            waiting += 1
    print(f"shared {len(shared)} | kept ours {mine} | linked {linked} (freed {freed / 1e9:.0f} GB) | waiting for the other job {waiting}")


if __name__ == "__main__":
    main()
