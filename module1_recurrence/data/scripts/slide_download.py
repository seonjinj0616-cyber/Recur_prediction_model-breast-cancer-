"""manifest.csv의 진단 슬라이드(환자당 첫 번째, 보통 DX1)를 받는다. 처음에는 300명부터.

디스크가 부족하면 feature 추출 뒤 svs를 지운다. svs나 feature(.h5)가 이미 있는 환자는 건너뛰므로
"받기 → 추출 → 삭제"를 묶음 단위로 반복할 수 있고, 중간에 끊겨도 다시 실행하면 된다.
{root}/qc/slide_qc.csv가 있으면 IHC로 판정된 슬라이드는 건너뛴다.
실행: python -m module1_recurrence.data.scripts.slide_download [--batch N] [--root DIR]
결과: {root}/slides/{file_name}.svs
"""
import argparse
import shutil
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from module1_recurrence.data.scripts.gdc_fetch import API, RAW

RETRIES = 3


def fetch(file_id, path):
    if path.exists():
        return "skip"
    tmp = path.with_suffix(".part")
    for attempt in range(RETRIES):
        try:
            # 슬라이드는 수 GB라서 메모리에 올리지 않고 바로 파일로 흘려 쓴다
            with urllib.request.urlopen(f"{API}/data/{file_id}") as r, open(tmp, "wb") as f:
                shutil.copyfileobj(r, f, length=1 << 20)
            tmp.rename(path)  # 다 받은 파일만 .svs 이름을 가진다
            return "ok"
        except OSError as e:  # GDC가 가끔 연결을 끊는다
            print(f"{path.name} retry {attempt + 1}: {e}")
            time.sleep(5 * (attempt + 1))
    return "failed"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", type=int, default=300, help="아직 처리하지 않은 환자 중 앞에서 N명")
    ap.add_argument("--root", type=Path, default=RAW.parent, help="manifest.csv가 있는 폴더 (서버: data/sj/tcga)")
    args = ap.parse_args()

    slide_dir, feat_dir = args.root / "slides", args.root / "features"
    slide_dir.mkdir(parents=True, exist_ok=True)

    m = pd.read_csv(args.root / "manifest.csv").dropna(subset=["slide_file_ids"])
    # manifest는 슬라이드를 파일 이름 순으로 묶어 두었으므로 첫 번째가 DX1
    m["file_id"] = m.slide_file_ids.str.split(";").str[0]
    m["path"] = [slide_dir / n for n in m.slide_file_names.str.split(";").str[0]]
    # slide_qc에서 IHC로 판정한 슬라이드(진단용 H&E가 아님)는 받지 않는다
    qc_path = args.root / "qc" / "slide_qc.csv"
    ihc = set(pd.read_csv(qc_path).query("stain == 'IHC_suspect'").file_id) if qc_path.exists() else set()
    m = m[~m.file_id.isin(ihc)]
    done = m.patient_id.map(lambda p: (feat_dir / f"{p}.h5").exists())
    todo = m[~done].head(args.batch)
    print(f"patients {len(m)} (IHC excluded {len(ihc)}) | features done {int(done.sum())} | this batch {len(todo)}")

    with ThreadPoolExecutor(max_workers=4) as ex:
        results = list(ex.map(fetch, todo.file_id, todo.path))

    gb = sum(f.stat().st_size for f in slide_dir.glob("*.svs")) / 1e9
    print({r: results.count(r) for r in set(results)}, f"| slides on disk {gb:.1f} GB")


if __name__ == "__main__":
    main()
