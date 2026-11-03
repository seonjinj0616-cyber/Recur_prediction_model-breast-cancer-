"""A3: 유사 RS가 METABRIC · SCAN-B에서 10년 원격 재발을 나누는지 확인한다.

TCGA와 같은 공식(pseudo_rs.normalize · score · assign_risk)을 각 데이터의 HR+/HER2− 전체 발현에 적용한다.
0–15 척도와 순위 3구간은 데이터마다 HR+/HER2− 전체 안에서 정한다 (마이크로어레이와 RNA-seq 척도가 달라서).
평가는 10년에서 자른 원격 재발: C-index(부트스트랩 95% CI), 표준편차당 Cox 위험비, 구간별 10년 KM 재발률, log-rank.
먼저 data.scripts.external_cohort, data.scripts.expression을 실행한다.
실행: python -m module1_recurrence.labels.scripts.pseudo_rs_prognosis
결과: data/{metabric,scanb}/pseudo_rs.csv, outputs/module1/a3_pseudo_rs_prognosis.json, a3_km_{metabric,scanb}.png
"""
import json

import matplotlib
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from lifelines.utils import concordance_index

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from module1_recurrence.data.scripts.expression import load  # noqa: E402
from module1_recurrence.data.scripts.external_fetch import DATA  # noqa: E402
from module1_recurrence.labels.scripts.pseudo_rs import assign_risk, normalize, score  # noqa: E402
from module1_recurrence.scripts.cohort import RECOMMENDED  # noqa: E402

OUT_DIR = DATA.parent / "outputs" / "module1"
TEN_YEARS = 3652.5
N_BOOT = 200
GROUPS = ["low", "intermediate", "high"]


def pseudo_rs(name):
    out = score(normalize(load(name)))
    out["rs_group"] = assign_risk(out.pseudo_rs)
    out.index.name = "patient_id"
    out.round(4).to_csv(DATA / name / "pseudo_rs.csv")
    return out


def subsets(name, df):
    """평가할 환자 범위. METABRIC은 재발 기록이 완전한 환자만 쓴다."""
    rec = df.group.isin(RECOMMENDED)
    if name == "metabric":
        df, rec = df[df.complete_rec], rec[df.complete_rec]
        endo_only = df.endocrine & ~df.chemo
    else:
        endo_only = df.treatment == "Endo"
    return {"all HR+/HER2-": df, "recommended": df[rec], "recommended + endocrine only": df[rec & endo_only]}


def truncate(df):
    event = (df.dr_event == 1) & (df.dr_days <= TEN_YEARS)
    return df.dr_days.clip(upper=TEN_YEARS), event.astype(int)


def evaluate(df, rng):
    t, e = truncate(df)
    rs = df.pseudo_rs
    cidx = concordance_index(t, rs * -1, e)  # 점수가 높을수록 빨리 재발
    boot = []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(df), len(df))
        boot.append(concordance_index(t.iloc[i], -rs.iloc[i], e.iloc[i]))

    cph = CoxPHFitter().fit(pd.DataFrame({"t": t, "e": e, "rs_sd": (rs - rs.mean()) / rs.std()}), "t", "e")
    s = cph.summary.loc["rs_sd"]

    km = {}
    for g in GROUPS:
        m = df.rs_group == g
        f = KaplanMeierFitter().fit(t[m], e[m])
        km[g] = {"n": int(m.sum()), "events": int(e[m].sum()),
                 "dr_10y_pct": round(100 * (1 - f.survival_function_at_times(TEN_YEARS - 1).iloc[0]), 1)}
    high = df.rs_group == "high"
    return {
        "n": len(df), "events_10y": int(e.sum()),
        "c_index": round(cidx, 3), "c_index_95ci": [round(float(x), 3) for x in np.percentile(boot, [2.5, 97.5])],
        "hr_per_sd": round(float(s["exp(coef)"]), 2),
        "hr_95ci": [round(float(s["exp(coef) lower 95%"]), 2), round(float(s["exp(coef) upper 95%"]), 2)], "hr_p": float(s.p),
        "km": km,
        "logrank_3group_p": float(multivariate_logrank_test(t, df.rs_group, e).p_value),
        "logrank_high_vs_rest_p": float(logrank_test(t[high], t[~high], e[high], e[~high]).p_value),
    }


def plot_km(name, subs, path):
    fig, axes = plt.subplots(1, len(subs), figsize=(5 * len(subs), 4), sharey=True)
    for ax, (label, df) in zip(axes, subs.items()):
        t, e = truncate(df)
        for g in GROUPS:
            m = df.rs_group == g
            KaplanMeierFitter(label=f"{g} (n={m.sum()})").fit(t[m] / 365.25, e[m]).plot_survival_function(
                ax=ax, ci_show=False)
        ax.set(title=f"{name} · {label}", xlabel="years", ylabel="distant recurrence-free")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    results = {}
    for name in ["metabric", "scanb"]:
        rs = pseudo_rs(name)
        cohort = pd.read_csv(DATA / name / "cohort.csv").dropna(subset=["dr_event", "dr_days"])
        df = cohort.join(rs[["pseudo_rs", "rs_group"]], on="patient_id")
        subs = subsets(name, df)
        results[name] = {label: evaluate(d.reset_index(drop=True), rng) for label, d in subs.items()}
        plot_km(name, subs, OUT_DIR / f"a3_km_{name}.png")

        print(f"\n{name}")
        for label, r in results[name].items():
            km = " / ".join(f"{g} {r['km'][g]['dr_10y_pct']}%" for g in GROUPS)
            print(f"  {label}: n={r['n']}, 10년 원격 재발 {r['events_10y']}건 | C-index {r['c_index']} "
                  f"{r['c_index_95ci']} | HR/SD {r['hr_per_sd']} {r['hr_95ci']} | 10년 재발률 {km} "
                  f"| log-rank 3구간 p={r['logrank_3group_p']:.1e}, high vs 나머지 p={r['logrank_high_vs_rest_p']:.1e}")

    (OUT_DIR / "a3_pseudo_rs_prognosis.json").write_text(json.dumps(results, ensure_ascii=False, indent=2),
                                                          encoding="utf-8")
    print("\nsaved", OUT_DIR)


if __name__ == "__main__":
    main()
