"""HistologyHSI-BC-Recurrence (TCIA, CC BY 4.0): H&E 슬라이드 47장과 임상 표만 받는다 (외부 원격 재발 검증).

TCIA는 1.2 TB(대부분 초분광 영상) 하나의 Aspera Faspex 패키지로만 준다. H&E 폴더(01_01_Histological_Images, 65 GB)와
조직 표시(01_02_Tissue_Annotations)만 골라 IBM aspera-cli(ascli)로 받는다 (2026-11-10, 실제로 받은 방법).
준비 (서버, 한 번):
  1. IBM Aspera Connect를 사용자 폴더에 설치 (~/.aspera/connect, 관리자 권한 불필요)
     https://d3gcli72yxqn2z.cloudfront.net/downloads/connect/latest/bin/ibm-aspera-connect_<버전>_linux_x86_64.tar.gz
  2. conda create -n aspera --override-channels -c conda-forge ruby=3.3 && conda activate aspera && gem install aspera-cli
     (Connect의 ascp 키는 비밀번호가 걸려 있어 ascp를 직접 부르면 막힘 → ascli가 인증을 처리함)
임상 표는 일반 HTTP로 받는다.
실행: python -m module1_recurrence.data.scripts.hsi_fetch --root data/sj/hsi
결과: {root}/01_01_Histological_Images/{id}.mrxs + {id}/, {root}/01_02_Tissue_Annotations/, {root}/clinical.xlsx
"""
import argparse
import json
import subprocess
import urllib.request
from pathlib import Path

LINK = ("https://faspex.cancerimagingarchive.net/aspera/faspex?context=eyJyZXNvdXJjZSI6InBhY2thZ2VzIiwidHlwZSI6ImV4dGVybmFsX2Rvd25sb2FkX3"
        "BhY2thZ2UiLCJpZCI6IjExNTAiLCJwYXNzY29kZSI6ImVhOGYxNDdhMmFlNzAyMjJiMmZlMGIyMjdkZjRhNGViMTNlMmI0YTAiLCJwYWNrYWdlX2lkIjoiMTE1"
        "MCIsImVtYWlsIjoiaGVscEBjYW5jZXJpbWFnaW5nYXJjaGl2ZS5uZXQifQ==")
FOLDERS = ["/01_01_Histological_Images", "/01_02_Tissue_Annotations"]
CLINICAL = "https://www.cancerimagingarchive.net/wp-content/uploads/HistologyHSI-BRCA-Recurrence-Clinical-Standardized.xlsx"
ASCLI = Path.home() / "miniconda3" / "envs" / "aspera" / "bin" / "ascli"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True, help="저장 폴더 (서버: data/sj/hsi)")
    args = ap.parse_args()
    args.root.mkdir(parents=True, exist_ok=True)
    xlsx = args.root / "clinical.xlsx"
    if not xlsx.exists():
        urllib.request.urlretrieve(CLINICAL, xlsx)
    ts = {"paths": [{"source": f} for f in FOLDERS], "resume_policy": "sparse_csum"}
    subprocess.run([str(ASCLI), "faspex5", "packages", "receive", f"--url={LINK}", "--sdk-folder=product:IBM Aspera Connect",
                    f"--to-folder={args.root.resolve()}", "--sources=@ts", f"--ts=@json:{json.dumps(ts)}", "--progress-bar=no"],
                   check=True)


if __name__ == "__main__":
    main()
