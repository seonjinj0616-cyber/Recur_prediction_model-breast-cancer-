"""특징 인코더 3개(ResNet50 · H-optimus-0 · Hibou-L)를 500장을 돌리기 전에 처리된 슬라이드로 비교한다 (B3 결정용).

슬라이드마다 RNA 조각 맞춤 타일과 무작위 조직 타일을 일부만 뽑아 한 번 읽고, 세 인코더에 같이 넣는다.
1. 속도: 인코더별 초당 타일 수 (GPU 계산만) + 타일 읽기 속도 → 처리된 슬라이드의 평균 조직 타일 수로 510장 시간 추정
2. 타일 미니 과제: CellViT 핵 구성(종양 · 염증 · 기질 비율, log 핵 수)을 ridge로 맞힘. 슬라이드 단위 5-fold, R²
3. 실제 과제 축소판: 환자별 맞춤 타일 특징 평균 → R 목록 · Oncotype 유전자 log2(TPM+1)를 ridge로 예측.
   환자 5-fold, 유전자별 Pearson r. 환자가 적어 흔들리므로 방향만 본다.
   신호가 없어도 교차검증 r은 0보다 낮게 나오므로(가짜 특징으로 시험 시 중앙값 −0.37), 같은 차원의 무작위 특징(random)을 기준으로 같이 계산한다.

H-optimus-0은 Hugging Face에서 이용 동의 후 `huggingface-cli login`이 필요하다 (사람이 직접).
Hibou-L은 모델 저장소의 코드를 실행한다 (trust_remote_code).

준비 (이 PC): python -m module1_recurrence.features.scripts.encoder_pilot --make_expr
  → data/tcga/expr_pilot.csv (환자 × 유전자)를 data/tcga/rna_purity.csv, outputs/module1/r_genes.csv와 함께 서버 {root}/로 복사
실행 (서버): python -m module1_recurrence.features.scripts.encoder_pilot --root data/sj/tcga [--n_slides 50] [--per_slide 300]
결과: {root}/pilot/ 아래 tiles.csv · {인코더}.npz · speed.json (같은 인코더 목록이면 캐시를 다시 씀), encoder_pilot.json
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from module1_recurrence.data.scripts.gdc_fetch import RAW
from module1_recurrence.features.scripts.tile_sets import select

ENCODERS = ["resnet50", "h_optimus_0", "hibou_l"]
TARGETS = ["tumor_ratio", "inflammatory_ratio", "connective_ratio", "log_cells"]


# ---------- 준비: 발현 표 (이 PC) ----------

def make_expr():
    from module1_recurrence.labels.scripts.pseudo_rs import CANCER_GENES
    genes = set(pd.read_csv(RAW.parents[2] / "outputs" / "module1" / "r_genes.csv").gene) | set(CANCER_GENES)
    m = pd.read_csv(RAW.parent / "manifest.csv")
    rows = {}
    for pid in m.patient_id:
        d = pd.read_csv(RAW.parent / "rnaseq" / f"{pid}.tsv", sep="\t", comment="#",
                        usecols=["gene_name", "tpm_unstranded"])
        rows[pid] = d[d.gene_name.isin(genes)].groupby("gene_name").tpm_unstranded.sum()
    expr = np.log2(pd.DataFrame(rows).T + 1)
    expr.index.name = "patient_id"
    expr.round(4).to_csv(RAW.parent / "expr_pilot.csv")
    print(expr.shape, "→", RAW.parent / "expr_pilot.csv")


# ---------- 인코더 ----------

def load_encoder(name, device):
    """(입력 NCHW float 0–1 → 특징) 함수와 정규화 평균 · 표준편차."""
    import torch
    if name == "resnet50":
        from module1_recurrence.features.scripts.extract import encoder
        return encoder(device), (0.485, 0.456, 0.406), (0.229, 0.224, 0.225)
    if name == "h_optimus_0":
        import timm
        net = timm.create_model("hf-hub:bioptimus/H-optimus-0", pretrained=True, init_values=1e-5,
                                dynamic_img_size=False).eval().to(device)
        return net, (0.707223, 0.578729, 0.703617), (0.211883, 0.230117, 0.177517)
    if name == "hibou_l":
        from transformers import AutoImageProcessor, AutoModel
        proc = AutoImageProcessor.from_pretrained("histai/hibou-L", trust_remote_code=True)
        model = AutoModel.from_pretrained("histai/hibou-L", trust_remote_code=True).eval().to(device)
        return (lambda x: model(pixel_values=x).pooler_output), tuple(proc.image_mean), tuple(proc.image_std)
    raise ValueError(name)


def read_tiles(svs, xy, tile, workers):
    """타일을 원본에서 읽어 224×224 uint8 (N, 224, 224, 3)."""
    from torch.utils.data import DataLoader
    from module1_recurrence.features.scripts.extract import Tiles

    class Raw(Tiles):
        def __getitem__(self, i):
            import openslide
            from PIL import Image
            if self.slide is None:
                self.slide = openslide.OpenSlide(self.svs)
            x, y = self.xy[i]
            img = self.slide.read_region((int(x), int(y)), 0, (self.tile, self.tile)).convert("RGB")
            return np.asarray(img.resize((224, 224), Image.BILINEAR))

    loader = DataLoader(Raw(svs, xy, tile), batch_size=64, num_workers=workers)
    return np.concatenate([b.numpy() for b in loader])


def embed(enc, imgs, device, batch_size):
    """특징 (N, D) float32와 GPU 계산 초."""
    import torch
    net, mean, std = enc
    mean = torch.tensor(mean, device=device).view(1, 3, 1, 1)
    std = torch.tensor(std, device=device).view(1, 3, 1, 1)
    out, sec = [], 0.0
    with torch.no_grad():
        for i in range(0, len(imgs), batch_size):
            x = torch.from_numpy(imgs[i:i + batch_size]).to(device).permute(0, 3, 1, 2).float() / 255
            x = (x - mean) / std
            torch.cuda.synchronize() if device.startswith("cuda") else None
            t0 = time.time()
            with torch.autocast(device_type=device.split(":")[0], dtype=torch.float16, enabled=device != "cpu"):
                f = net(x)
            torch.cuda.synchronize() if device.startswith("cuda") else None
            sec += time.time() - t0
            out.append(f.float().cpu().numpy())
    return np.concatenate(out), sec


# ---------- 타일 뽑기 · 특징 ----------

def sample_tiles(csv, purity, per_slide, rng):
    """맞춤 타일 최대 per_slide개 + 그 밖의 조직 타일 최대 per_slide개."""
    t = pd.read_csv(csv)
    t = t[t.tissue_tile].reset_index(drop=True)
    matched = select(t.tumor_ratio, t.n_cells, "matched", purity=purity)
    a = np.flatnonzero(matched)
    b = np.flatnonzero(~matched)
    pick = np.r_[rng.permutation(a)[:per_slide], rng.permutation(b)[:per_slide]]
    s = t.loc[pick].reset_index(drop=True)
    s["matched"] = matched[pick]
    return s, len(t)


def collect(root, slides, purity, args):
    """모든 인코더 특징을 만들거나 캐시에서 읽는다 → (타일 표, {인코더: 특징}, 속도)."""
    out = root / "pilot"
    out.mkdir(exist_ok=True)
    meta_path, speed_path = out / "tiles.csv", out / "speed.json"
    todo = args.encoders
    if speed_path.exists() and set(todo) <= set(json.loads(speed_path.read_text())["encoders"]):
        return (pd.read_csv(meta_path), {e: np.load(out / f"{e}.npz")["x"] for e in todo},
                json.loads(speed_path.read_text()))

    rng = np.random.default_rng(0)
    encs = {e: load_encoder(e, args.device) for e in todo}
    tables, feats, gpu_sec, read_sec, n_tissue = [], {e: [] for e in todo}, dict.fromkeys(todo, 0.0), 0.0, []
    for slide, pid, svs in slides:
        s, n = sample_tiles(root / "tiles" / f"{slide}.csv", purity.get(pid), args.per_slide, rng)
        n_tissue.append(n)
        t0 = time.time()
        imgs = read_tiles(svs, s[["x", "y"]].to_numpy(), s.tile.iloc[0], args.workers)
        read_sec += time.time() - t0
        for e in todo:
            f, sec = embed(encs[e], imgs, args.device, args.batch_size)
            feats[e].append(f)
            gpu_sec[e] += sec
        s.insert(0, "slide", slide)
        s.insert(1, "patient_id", pid)
        tables.append(s)
        print(slide, len(s), "tiles", {e: round(gpu_sec[e], 1) for e in todo})

    meta = pd.concat(tables, ignore_index=True)
    meta.to_csv(meta_path, index=False)
    for e in todo:
        np.savez(out / f"{e}.npz", x=np.concatenate(feats[e]).astype(np.float16))
    n = len(meta)
    speed = {"tiles": n, "read_tiles_per_s": round(n / read_sec, 1),
             "mean_tissue_tiles_per_slide": round(float(np.mean(n_tissue))),
             "encoders": {e: {"tiles_per_s": round(n / gpu_sec[e], 1),
                              # 510장 × 평균 조직 타일, 읽기와 계산 중 느린 쪽이 전체 속도를 정함
                              "hours_510_slides": round(510 * np.mean(n_tissue) / min(n / gpu_sec[e], n / read_sec) / 3600, 1)}
                          for e in todo}}
    speed_path.write_text(json.dumps(speed, indent=2))
    return meta, {e: np.load(out / f"{e}.npz")["x"] for e in args.encoders}, speed


# ---------- 평가 ----------

def tile_targets(meta):
    n = meta.n_cells.to_numpy(dtype=float)
    return pd.DataFrame({"tumor_ratio": meta.tumor_ratio, "inflammatory_ratio": meta.n_inflammatory / n,
                         "connective_ratio": meta.n_connective / n, "log_cells": np.log(n)})


def probe_tiles(x, meta):
    """슬라이드 단위 5-fold ridge → 과제별 R²."""
    from sklearn.linear_model import RidgeCV
    from sklearn.metrics import r2_score
    from sklearn.model_selection import GroupKFold, cross_val_predict
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    y = tile_targets(meta)
    model = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-1, 5, 13)))
    pred = cross_val_predict(model, x.astype(np.float32), y.to_numpy(), groups=meta.slide, cv=GroupKFold(5))
    return {c: round(r2_score(y[c], pred[:, i]), 3) for i, c in enumerate(TARGETS)}


def probe_expression(x, meta, expr, r_genes, onco):
    """환자별 맞춤 타일 특징 평균 → 유전자 log2(TPM+1), 환자 5-fold ridge → 유전자별 Pearson r 요약."""
    from sklearn.linear_model import RidgeCV
    from sklearn.model_selection import KFold, cross_val_predict
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    m = meta.matched.to_numpy()
    pat = pd.DataFrame(x[m].astype(np.float32)).groupby(meta.patient_id[m].to_numpy()).mean()
    pat = pat.loc[pat.index.isin(expr.index)]
    y = expr.loc[pat.index]
    y = y.loc[:, y.std() > 0]
    model = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(0, 6, 13)))
    pred = cross_val_predict(model, pat.to_numpy(), y.to_numpy(), cv=KFold(5, shuffle=True, random_state=0))
    yc, pc = y.to_numpy() - y.to_numpy().mean(0), pred - pred.mean(0)
    r = pd.Series((yc * pc).sum(0) / np.sqrt((yc ** 2).sum(0) * (pc ** 2).sum(0) + 1e-12), index=y.columns)

    def summary(g):
        v = r[r.index.isin(g)]
        return {"genes": len(v), "median_r": round(float(v.median()), 3), "r_gt_0.3": int((v > 0.3).sum())}
    return {"patients": len(pat), "r_genes": summary(r_genes), "oncotype": summary(onco)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--make_expr", action="store_true", help="이 PC에서 data/tcga/expr_pilot.csv만 만든다")
    ap.add_argument("--root", type=Path, default=RAW.parent, help="slides/ · tiles/ · rna_purity.csv · expr_pilot.csv · r_genes.csv 폴더")
    ap.add_argument("--n_slides", type=int, default=50)
    ap.add_argument("--per_slide", type=int, default=300, help="슬라이드마다 맞춤 타일 · 그 밖의 조직 타일을 각각 최대 이만큼")
    ap.add_argument("--encoders", nargs="+", default=ENCODERS, choices=ENCODERS)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--batch_size", type=int, default=64)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    if args.make_expr:
        return make_expr()

    from module1_recurrence.labels.scripts.pseudo_rs import CANCER_GENES
    root = args.root
    purity = (pd.read_csv(root / "rna_purity.csv").set_index("patient_id").tumor_nuclei / 100).to_dict()
    expr = pd.read_csv(root / "expr_pilot.csv", index_col="patient_id")
    svs = {p.name.split(".")[0]: p for p in (root / "slides").glob("*.svs")}
    slides = [(c.stem, c.stem[:12], svs[c.stem]) for c in sorted((root / "tiles").glob("TCGA-*.csv"))
              if c.stem in svs and c.stem[:12] in purity][:args.n_slides]
    print(len(slides), "slides")

    meta, feats, speed = collect(root, slides, purity, args)
    r_genes = set(pd.read_csv(root / "r_genes.csv").gene)
    result = {"slides": meta.slide.nunique(), "tiles": len(meta), "speed": speed, "encoders": {}}
    feats["random"] = np.random.default_rng(0).normal(size=(len(meta), 1024)).astype(np.float16)  # 기준선
    for e, x in feats.items():
        result["encoders"][e] = {"dim": int(x.shape[1]), "tile_r2": probe_tiles(x, meta),
                                 "expression": probe_expression(x, meta, expr, r_genes, set(CANCER_GENES))}
        print(e, json.dumps(result["encoders"][e], ensure_ascii=False))
    (root / "pilot" / "encoder_pilot.json").write_text(json.dumps(result, indent=2, ensure_ascii=False))
    print(json.dumps(speed, indent=2))


if __name__ == "__main__":
    main()
