"""H&E → 원격 재발 평가 (로컬): H&E 점수(CellViT 없음, 전체 조직)가 10년 원격 재발을 나누는지 본다. 원격 재발만 쓴다.

입력: 서버의 c2_validation_scores.csv와 HistologyHSI-BC 임상 표(clinical.xlsx)를 data/hsi/에 둔다.
- 주 평가 (모듈 대상 HR+/HER2−만)
  · TCGA HR+/HER2−: 기관 단위 교차검증 점수, 10년 원격 전이 (일). 측정 RNA 점수 C-index를 상한선으로 같이 계산.
  · HistologyHSI-BC HR+/HER2−: 원격 재발 Relapse, 시간 DFS (개월). 재발군 · 무재발군을 골라 모은 환자-대조 설계라 순위 지표만 해석.
  · 합산: 두 데이터는 시간 단위 · 추적이 달라 같은 데이터 안의 환자 쌍만 비교하는 층화 C-index
    (데이터별 일치 쌍 수 합 / 비교 가능 쌍 수 합). 부트스트랩은 데이터마다 따로 다시 뽑는다.
- 보조: TCGA 검사 권고 코호트(N0 · N1 폐경 후), HistologyHSI-BC 전체 47명(삼중음성 · HER2 양성 포함).
비교 기준: NPI (HistologyHSI-BC는 크기 · 등급 · 림프절로 계산, TCGA는 등급이 없어 계산 안 함).
실행: python -m module1_recurrence.models.scripts.c2_hsi_eval
결과: outputs/module1/c2_validation.json
"""
import json

import numpy as np
import pandas as pd

from module1_recurrence.data.scripts.external_fetch import DATA

OUT_DIR = DATA.parent / "outputs" / "module1"
HSI = DATA / "hsi"
N_BOOT = 1000


def pairs(t, e, s):
    """일치 쌍 수와 비교 가능 쌍 수 (점수가 높을수록 빨리 재발이 정답, 동점은 0.5)."""
    t, e, s = map(np.asarray, (t, e, s))
    conc = comp = 0.0
    for i in np.flatnonzero(e == 1):
        later = (t > t[i]) | ((t == t[i]) & (e == 0))  # 같은 시점의 비사건 환자도 비교 (lifelines와 같음)
        comp += later.sum()
        conc += (s[i] > s[later]).sum() + 0.5 * (s[i] == s[later]).sum()
    return conc, comp


def stratified(groups):
    c = [pairs(*g) for g in groups]
    return sum(a for a, _ in c) / sum(b for _, b in c)


def boot_ci(groups, rng):
    vals = []
    for _ in range(N_BOOT):
        g2 = []
        for t, e, s in groups:
            i = rng.integers(0, len(t), len(t))
            g2.append((np.asarray(t)[i], np.asarray(e)[i], np.asarray(s)[i]))
        vals.append(stratified(g2))
    return [round(float(v), 3) for v in np.percentile(vals, [2.5, 97.5])]


def summarize(groups, rng):
    return {"n": int(sum(len(g[0]) for g in groups)), "events": int(sum(np.sum(g[1]) for g in groups)),
            "c_index": round(stratified(groups), 3), "95ci": boot_ci(groups, rng)}


def main():
    rng = np.random.default_rng(0)
    s = pd.read_csv(HSI / "c2_validation_scores.csv", index_col="id", dtype={"id": str})
    clin = pd.read_excel(HSI / "clinical.xlsx", sheet_name="Clinical Data")
    clin.index = clin["Case ID"].astype(str)
    hsi = s[s.cohort == "HistologyHSI-BC"].drop(columns=["time", "event"], errors="ignore").join(clin, how="inner")
    hsi["time"], hsi["event"] = hsi.DFS, hsi.Relapse.astype(int)
    hsi["npi"] = 0.2 * hsi.Tumor_diameter / 10 + hsi.Tumor_histologic_grade + hsi.N.map({0: 1, 1: 2, 2: 3, 3: 3})
    hsi_hr = hsi[((hsi.ER == 1) | (hsi.PR == 1)) & (hsi.HER2 == 0)]
    tcga = s[s.cohort == "TCGA HR+/HER2- (CV)"]
    tcga_rec = tcga[tcga.recommended.astype(str) == "True"]

    g = lambda d, col="score": (d.time.to_numpy(), d.event.to_numpy(), d[col].to_numpy())  # noqa: E731
    res = {"tcga_cv_r_with_rna_score": round(float(np.corrcoef(tcga.score, tcga.measured)[0, 1]), 3), "primary (HR+/HER2-)": {
        "TCGA (CV)": {**summarize([g(tcga)], rng), "measured_rna_score": summarize([g(tcga, "measured")], rng)},
        "HistologyHSI-BC": {**summarize([g(hsi_hr)], rng), "npi": summarize([g(hsi_hr, "npi")], rng)},
        "pooled (stratified)": summarize([g(tcga), g(hsi_hr)], rng)}, "secondary": {
        "TCGA recommended (CV)": {**summarize([g(tcga_rec)], rng), "measured_rna_score": summarize([g(tcga_rec, "measured")], rng)},
        "HistologyHSI-BC all 47": {**summarize([g(hsi)], rng), "npi": summarize([g(hsi, "npi")], rng)}}}
    (OUT_DIR / "c2_validation.json").write_text(json.dumps(res, indent=2, ensure_ascii=False))
    print(json.dumps(res, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
