"""ABMIL (방법 3): 타일 특징 → gated attention pooling (Ilse et al. 2018) → C2 재발 점수 회귀.

평균 대신 타일마다 중요도(attention)를 학습해 가중 평균한다. 입력은 환자의 맞춤 타일 (tile_sets.select, p_genes와 같음).
정답은 c2_he_direct가 만든 측정 RNA C2 점수 (c2_he/c2_he_scores.csv의 truth). 기관(TSS) 단위 5-fold 교차검증 예측을 만들어
평균 방식(direct)과 같은 지표(정답 점수와의 r, 10년 PFI C-index)로 비교한다.
학습: 정답을 표준화, MSE, AdamW, EPOCHS번, 가방당 최대 MAX_TILES개 무작위 타일. 시드 SEEDS개 평균.
실행: python -m module1_recurrence.models.scripts.abmil --root data/sj/tcga --gpu 0
결과: {root}/c2_he/abmil_scores.csv, abmil.json
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from lifelines.utils import concordance_index
from sklearn.model_selection import GroupKFold
from torch import nn

from module1_recurrence.features.scripts.tile_sets import select

EPOCHS = 30
MAX_TILES = 2000
SEEDS = 3
HIDDEN = 256


class ABMIL(nn.Module):
    def __init__(self, dim=1536, hidden=HIDDEN, dropout=0.25):
        super().__init__()
        self.embed = nn.Sequential(nn.Linear(dim, hidden), nn.GELU(), nn.Dropout(dropout))
        self.att_v = nn.Sequential(nn.Linear(hidden, hidden // 2), nn.Tanh())
        self.att_u = nn.Sequential(nn.Linear(hidden, hidden // 2), nn.Sigmoid())
        self.att_w = nn.Linear(hidden // 2, 1)
        self.head = nn.Linear(hidden, 1)

    def forward(self, x):  # x: (타일, dim)
        h = self.embed(x)
        a = torch.softmax(self.att_w(self.att_v(h) * self.att_u(h)).squeeze(-1), dim=0)
        return self.head((a[:, None] * h).sum(0)).squeeze(-1), a


def load_bags(root, pids):
    import h5py
    purity = (pd.read_csv(root / "rna_purity.csv").set_index("patient_id").tumor_nuclei / 100).to_dict()
    bags = {}
    for pid in pids:
        f = next((root / "features").glob(f"{pid}-*DX1.h5"))
        with h5py.File(f, "r") as h:
            keep = select(h["tumor_ratio"][:], h["n_cells"][:], "matched", purity=purity.get(pid))
            bags[pid] = torch.from_numpy(h["features"][:][keep].astype(np.float32))
    return bags


def train_predict(bags, y, tr, te, device, seed):
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    model = ABMIL().to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=1e-4, weight_decay=1e-2)
    for _ in range(EPOCHS):
        model.train()
        for i in rng.permutation(tr):
            x = bags[i]
            if len(x) > MAX_TILES:
                x = x[rng.choice(len(x), MAX_TILES, replace=False)]
            loss = (model(x.to(device))[0] - y[i]) ** 2
            opt.zero_grad()
            loss.backward()
            opt.step()
    model.eval()
    with torch.no_grad():
        return np.array([model(bags[i].to(device))[0].item() for i in te])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--gpu", default="0")
    args = ap.parse_args()
    device = f"cuda:{args.gpu}" if torch.cuda.is_available() else "cpu"
    s = pd.read_csv(args.root / "c2_he" / "c2_he_scores.csv", index_col=0)
    pids = s.index.to_list()
    bags = load_bags(args.root, pids)
    bag_list = [bags[p] for p in pids]
    y = ((s.truth - s.truth.mean()) / s.truth.std()).to_numpy(dtype=np.float32)
    groups = s.index.str[5:7].to_numpy()
    pred = np.zeros(len(pids))
    for fold, (tr, te) in enumerate(GroupKFold(5).split(pids, groups=groups)):
        p = np.mean([train_predict(bag_list, torch.tensor(y), tr, te, device, seed) for seed in range(SEEDS)], axis=0)
        pred[te] = p
        print(f"fold {fold}: r {np.corrcoef(p, y[te])[0, 1]:.3f}", flush=True)
    s["abmil"] = pred
    s[["truth", "abmil", "T", "E"]].to_csv(args.root / "c2_he" / "abmil_scores.csv")
    t, e = s["T"].to_numpy(), s["E"].to_numpy().astype(int)
    res = {k: {"r_with_truth": round(float(np.corrcoef(s[k], s.truth)[0, 1]), 3),
               "pfi_c_index": round(concordance_index(t, -s[k].to_numpy(), e), 3)} for k in ["direct", "gene", "abmil"]}
    (args.root / "c2_he" / "abmil.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
