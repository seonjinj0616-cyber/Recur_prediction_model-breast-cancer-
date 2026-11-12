"""TCGA HR+/HER2−: 원격 전이 기록(new tumor event)이 다른 재발 관련 기록과 맞는지 대조하고, 기록 누락이 C-index에 주는 영향을 본다.

원격 재발 평가는 임상 new tumor event 중 Distant Metastasis만 쓴다 (data/tcga/distant.csv, manifest.load_distant).
이 기록은 병원마다 추가 양식으로 적어 빠진 경우가 있어, 아래와 대조한다.
- TCGA-CDR PFI(국소 재발 · 원격 전이 · 새 원발암 · 진행 · 암 사망)와 DSS(유방암 사망)
- new tumor event의 종류가 국소 재발 · 새 원발암인지, 비어 있는지
- 추적 시 종양 상태(WITH TUMOR)
기록 누락 의심 = 원격 전이 기록 없이 (종류가 비어 있거나 국소 · 새 원발암이 아닌 PFI 사건) 또는 (유방암 사망).
민감도 분석: 의심 사례를 원격 재발로 넣고(시간은 PFI 또는 DSS 일수) H&E · 측정 RNA 점수 C-index를 다시 계산.
의심 사례는 확정할 수 없으므로 본 결과가 아니라 민감도 분석으로만 쓴다.
먼저 manifest(distant.csv)와 c2_hsi_scores · c2_hsi_eval(data/hsi/c2_validation_scores.csv)을 실행한다.
실행: python -m module1_recurrence.data.scripts.tcga_distant_audit
결과: outputs/module1/tcga_distant_audit.json
"""
import json

import numpy as np
import pandas as pd
from lifelines.utils import concordance_index

from module1_recurrence.data.scripts.external_fetch import DATA
from module1_recurrence.data.scripts.gdc_fetch import RAW
from module1_recurrence.data.scripts.manifest import NTE_FILES

OUT = DATA.parent / "outputs" / "module1" / "tcga_distant_audit.json"
FOLLOW_UP = ["nationwidechildrens.org_clinical_follow_up_v1.5_brca.txt", "nationwidechildrens.org_clinical_follow_up_v2.1_brca.txt",
             "nationwidechildrens.org_clinical_follow_up_v4.0_brca.txt"]
TEN_YEARS = 3652.5


def main():
    m = pd.read_csv(RAW.parent / "manifest.csv").set_index("patient_id")
    ids = m.index
    d = pd.read_csv(RAW.parent / "distant.csv", index_col=0).reindex(ids)
    cdr = pd.read_excel(RAW / "TCGA-CDR-SupplementalTableS1.xlsx", sheet_name="TCGA-CDR")
    cdr = cdr[cdr.type == "BRCA"].set_index("bcr_patient_barcode").reindex(ids)

    nte = pd.concat([pd.read_csv(RAW / f, sep="\t", skiprows=[1, 2])[["bcr_patient_barcode", "new_tumor_event_type"]] for f in NTE_FILES])
    nte = nte[nte.bcr_patient_barcode.isin(ids)]
    types = nte.new_tumor_event_type.astype(str)
    local_or_new = set(nte[types.str.contains("Locoregional|New Primary")].bcr_patient_barcode)
    type_missing = set(nte[types.str.contains("Not Available|Unknown")].bcr_patient_barcode)
    fu = pd.concat([pd.read_csv(RAW / f, sep="\t", skiprows=[1, 2])[["bcr_patient_barcode", "tumor_status"]] for f in FOLLOW_UP])
    with_tumor = set(fu[fu.tumor_status == "WITH TUMOR"].bcr_patient_barcode) & set(ids)

    dm, pfi, dss = d.dm_event == 1, cdr.PFI == 1, cdr.DSS == 1
    pfi_no_dm = pfi & ~dm
    cand = (pfi_no_dm & ~ids.isin(local_or_new)) | (dss & ~dm)
    res = {"patients": len(ids), "distant_recorded": int(dm.sum()), "pfi_events": int(pfi.sum()), "dss_events": int(dss.sum()),
           "pfi_without_distant": int(pfi_no_dm.sum()),
           "pfi_without_distant_local_or_new_primary": len(local_or_new & set(ids[pfi_no_dm])),
           "pfi_without_distant_type_missing": len(type_missing & set(ids[pfi_no_dm])),
           "dss_death_without_distant": int((dss & ~dm).sum()),
           "follow_up_with_tumor_without_distant": len(with_tumor - set(ids[dm])),
           "suspected_unrecorded_distant": int(cand.sum()),
           "sites_with_followup_over_5y": f"{int((d.dm_days > 1826).groupby(m.tss).any().sum())} of {m.tss.nunique()}",
           "followup_median_years": round(float(d.dm_days.median() / 365.25), 1)}

    s = pd.read_csv(DATA / "hsi" / "c2_validation_scores.csv", index_col="id")
    s = s[s.cohort == "TCGA HR+/HER2- (CV)"]
    t2 = np.where(cand, np.where(dss & ~dm, cdr["DSS.time"], cdr["PFI.time"]), d.dm_days)
    o = pd.DataFrame({"t2": np.minimum(t2, TEN_YEARS), "e2": (dm | cand).astype(int).to_numpy() * (t2 <= TEN_YEARS)}, index=ids)
    o = o.join(s[["score", "measured", "time", "event"]], how="inner")
    res["sensitivity"] = {lab: {"events_10y": int(e.sum()),
                                "he_c_index": round(concordance_index(t, -o.score, e), 3),
                                "rna_c_index": round(concordance_index(t, -o.measured, e), 3)}
                          for lab, t, e in [("recorded distant", o.time, o.event), ("+ suspected unrecorded", o.t2, o.e2)]}
    OUT.write_text(json.dumps(res, indent=2, ensure_ascii=False))
    print(json.dumps(res, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
