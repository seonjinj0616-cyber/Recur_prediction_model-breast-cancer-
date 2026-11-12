"""외부 검증 (로컬): 학습에 쓰지 않은 환자의 H&E 점수(CellViT 없음, 전체 조직)가 재발을 나누는지 본다.

입력: 서버의 c2_validation_scores.csv와 HistologyHSI-BC 임상 표(clinical.xlsx)를 data/hsi/에 둔다.
- HistologyHSI-BC 47명: 재발 = 원격 재발 Relapse, 시간 = DFS (개월). 재발군 22 · 무재발군 25를 골라 모은 환자-대조 설계라
  재발률은 실제보다 높고 순위 지표만 해석한다. HR+/HER2− 33명 따로.
- TCGA 다른 아형 (학습 · 유전자 선택에 안 씀, 진단 시 M1 제외): 재발 = 10년 원격 전이 (일, 임상 new tumor event).
  PFI(국소 재발 · 진행 포함) C-index도 참고로 계산. 대부분 HER2 양성 · 삼중음성이라 HR+/HER2− 환자만 따로도 봄.
  측정 RNA 점수의 C-index를 상한선으로 같이 계산.
합친 통계: 두 데이터는 재발 정의 · 시간 단위가 달라, 같은 데이터 안의 환자 쌍만 비교하는 층화 C-index
(데이터별 일치 쌍 수 합 / 비교 가능 쌍 수 합). 부트스트랩은 데이터마다 따로 다시 뽑는다.
비교 기준: NPI (HistologyHSI-BC는 크기 · 등급 · 림프절로 계산, TCGA는 등급이 없어 계산 안 함).
실행: python -m module1_recurrence.models.scripts.c2_hsi_eval
결과: outputs/module1/c2_validation.json
"""
import json

import numpy as np
import pandas as pd
from lifelines.utils import concordance_index

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


def main():
    rng = np.random.default_rng(0)
    s = pd.read_csv(HSI / "c2_validation_scores.csv", index_col="id", dtype={"id": str})
    clin = pd.read_excel(HSI / "clinical.xlsx", sheet_name="Clinical Data")
    clin.index = clin["Case ID"].astype(str)
    hsi = s[s.cohort == "HistologyHSI-BC"].join(clin, how="inner")
    hsi["time"], hsi["event"] = hsi.DFS, hsi.Relapse.astype(int)
    hsi["npi"] = 0.2 * hsi.Tumor_diameter / 10 + hsi.Tumor_histologic_grade + hsi.N.map({0: 1, 1: 2, 2: 3, 3: 3})
    hr = ((hsi.ER == 1) | (hsi.PR == 1)) & (hsi.HER2 == 0)
    other = s[s.cohort == "TCGA other subtypes"]
    other_hr = other[other.group != "not_HR+HER2-"]

    sets = {"HistologyHSI-BC (all)": hsi, "HistologyHSI-BC (HR+/HER2-)": hsi[hr],
            "TCGA other subtypes (all)": other, "TCGA other subtypes (HR+/HER2-)": other_hr}
    res = {}
    for k, d in sets.items():
        if len(d) == 0 or d.event.sum() == 0:
            res[k] = {"n": len(d), "events": int(d.event.sum()) if len(d) else 0}
            continue
        g = [(d.time, d.event, d.score)]
        res[k] = {"n": len(d), "events": int(d.event.sum()), "c_index": round(stratified(g), 3), "95ci": boot_ci(g, rng)}
        if "pfi_event" in d and d.pfi_event.notna().all():
            gp = [(d.pfi_time, d.pfi_event, d.score)]
            res[k]["pfi_c_index"] = {"c_index": round(stratified(gp), 3), "95ci": boot_ci(gp, rng), "events": int(d.pfi_event.sum())}
        for ref in ["npi", "measured"]:
            if ref in d and d[ref].notna().all():
                gr = [(d.time, d.event, d[ref])]
                res[k][ref] = {"c_index": round(stratified(gr), 3), "95ci": boot_ci(gr, rng)}
        if "measured" in d and d.measured.notna().all():
            res[k]["r_he_vs_measured"] = round(float(np.corrcoef(d.score, d.measured)[0, 1]), 3)
    for k, (a, b) in {"pooled HistologyHSI-BC + TCGA other (all)": (hsi, other),
                      "pooled HistologyHSI-BC + TCGA other (HR+/HER2-)": (hsi[hr], other_hr)}.items():
        groups = [(d.time, d.event, d.score) for d in (a, b) if len(d) and d.event.sum() > 0]
        res[k] = {"n": len(a) + len(b), "events": int(a.event.sum() + b.event.sum()),
                  "stratified_c_index": round(stratified(groups), 3), "95ci": boot_ci(groups, rng)}
    (OUT_DIR / "c2_validation.json").write_text(json.dumps(res, indent=2, ensure_ascii=False))
    print(json.dumps(res, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
