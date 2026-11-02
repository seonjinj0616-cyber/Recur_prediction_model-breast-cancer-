"""METABRIC · SCAN-B 임상 + 원격 재발을 환자 1행 코호트 표로 만든다.

먼저 external_fetch를 실행한다.
METABRIC 원격 재발: Rueda et al. 2019 Table S5 (DR, TDR 일수). SCAN-B 원격 재발: DRFi (일수).
실행: python -m module1_recurrence.data.external_cohort
결과: data/metabric/cohort.csv (HR+/HER2− 1,413명), data/scanb/cohort.csv (HR+/HER2− 추적 코호트)
"""
import pandas as pd

from module1_recurrence.cohort import RECOMMENDED, standardize_metabric, standardize_scanb
from module1_recurrence.data.external_fetch import DATA, METABRIC, SCANB

TEN_YEARS = 3652.5


def build_metabric():
    clin = pd.read_csv(DATA / "metabric" / "hr_her2neg.csv")
    rueda = pd.read_csv(METABRIC / "rueda2019_TableS5.txt", sep="\t").rename(columns={"METABRIC.ID": "patient_id"})

    out = standardize_metabric(clin).assign(
        age=clin.AGE_AT_DIAGNOSIS, tumor_size_mm=clin.TUMOR_SIZE, grade=clin.GRADE, npi=clin.NPI,
        chemo=clin.CHEMOTHERAPY == "YES", endocrine=clin.HORMONE_THERAPY == "YES",
        expr_id=clin.patientId)
    out = out.merge(rueda[["patient_id", "DR", "TDR", "LR", "TLR", "DeathBreast", "Death", "T",
                           "Complete.Rec.History"]], on="patient_id", how="left")
    out = out.rename(columns={"DR": "dr_event", "TDR": "dr_days", "Complete.Rec.History": "complete_rec"})
    out["complete_rec"] = out.complete_rec == "YES"
    return out


def build_scanb():
    clin = pd.read_excel(SCANB / "scanb_clinical_2023-01-13.xlsx", sheet_name="SCANB.9206")
    clin = clin[clin["Follow.up.cohort"]].reset_index(drop=True)  # 환자당 1샘플

    out = standardize_scanb(clin).assign(
        age_range=clin[clin.columns[clin.columns.str.startswith("Age")][0]], tumor_size_mm=clin["Size.mm"],
        grade=clin.NHG, treatment=clin.TreatGroup,
        chemo=clin.TreatGroup.str.contains("Chemo", na=False), endocrine=clin.TreatGroup.str.contains("Endo", na=False),
        expr_id=clin["GEX.assay"],
        dr_event=clin.DRFi_event, dr_days=clin.DRFi_days)
    return out[out.group != "not_HR+HER2-"].reset_index(drop=True)


def summarize(df, name):
    df = df.dropna(subset=["dr_event"])
    g = df.groupby("group")
    s = pd.DataFrame({"n": g.size(), "dr": g.dr_event.sum().astype(int),
                      "dr_10y": g.apply(lambda x: int(((x.dr_event == 1) & (x.dr_days <= TEN_YEARS)).sum())),
                      "median_fu_y": (g.dr_days.median() / 365.25).round(1)})
    print(f"\n{name} (원격 재발 정보 있음 {len(df)}명)")
    print(s)
    rec = df[df.group.isin(RECOMMENDED)]
    print(f"검사 권고 {len(rec)}명, 원격 재발 {int(rec.dr_event.sum())}건"
          f" | 내분비요법만 {(rec.endocrine & ~rec.chemo).sum()}명, {int(rec[rec.endocrine & ~rec.chemo].dr_event.sum())}건")
    return rec


def main():
    mb = build_metabric()
    mb.to_csv(DATA / "metabric" / "cohort.csv", index=False)
    summarize(mb, "METABRIC HR+/HER2−")
    rec = summarize(mb[mb.complete_rec], "METABRIC 재발 기록 완전")
    print("치료 (검사 권고 + 기록 완전):", rec.groupby(["chemo", "endocrine"]).size().to_dict())

    sb = build_scanb()
    sb.to_csv(DATA / "scanb" / "cohort.csv", index=False)
    summarize(sb, "SCAN-B HR+/HER2−")
    print("\nsaved", DATA / "metabric" / "cohort.csv", "|", DATA / "scanb" / "cohort.csv")


if __name__ == "__main__":
    main()
