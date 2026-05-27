#!/usr/bin/env python3
"""
annotate_sex_check.py
Merges PLINK --check-sex output with collected sex from sex_info.
Flags samples where inferred sex disagrees with collected sex.

Changes from previous version (H3AGWAS alignment):
  - F-statistic thresholds are now configurable via --f_lo_male / --f_hi_female
    rather than hardcoded 0.2 / 0.8, matching H3AGWAS params f_lo_male / f_hi_female.
  - Ambiguous range is anything strictly between the two thresholds.

PLINK sexcheck columns:
  FID IID PEDSEX SNPSEX STATUS F
  SNPSEX: 1=male, 2=female, 0=unknown
  STATUS: OK / PROBLEM
"""
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--sexcheck",    required=True)
parser.add_argument("--sex_info",    required=True)   # sampleid, sex (0=F, 1=M)
parser.add_argument("--out_annot",   required=True)
parser.add_argument("--out_discord", required=True)
# H3AGWAS-derived thresholds (defaults match PLINK built-ins)
parser.add_argument("--f_lo_male",   type=float, default=0.8,
                    help="F-statistic lower bound for calling Male (default: 0.8)")
parser.add_argument("--f_hi_female", type=float, default=0.2,
                    help="F-statistic upper bound for calling Female (default: 0.2)")
args = parser.parse_args()

f_lo_male   = args.f_lo_male
f_hi_female = args.f_hi_female

print(f"F-statistic thresholds: Male > {f_lo_male}, Female < {f_hi_female}")

sexcheck = pd.read_csv(args.sexcheck, sep=r"\s+")

# IIDs from PLINK are read as int64 when sample IDs are numeric.
# Cast both sides to str so the merge key types match.
sexcheck["IID"] = sexcheck["IID"].astype(str).str.strip()
sexcheck["FID"] = sexcheck["FID"].astype(str).str.strip()

sex_info = pd.read_csv(args.sex_info, sep="\t", dtype=str)
sex_info.columns = sex_info.columns.str.strip().str.lower()
sex_info = sex_info.rename(columns={"sampleid": "IID"})
sex_info["IID"] = sex_info["IID"].astype(str).str.strip()

# Map collected sex: 0→Female, 1→Male
sex_info["COLLECTED_SEX"] = sex_info["sex"].map({"0": "Female", "1": "Male"}).fillna("Unknown")

# Map PLINK inferred sex: 1→Male, 2→Female, 0→Unknown
sexcheck["INFERRED_SEX"] = sexcheck["SNPSEX"].map({1: "Male", 2: "Female", 0: "Unknown"})

# F-statistic label using configurable thresholds (H3AGWAS alignment)
def f_label(f):
    if pd.isna(f):         return "Unknown"
    if f > f_lo_male:      return "Male"
    if f < f_hi_female:    return "Female"
    return "Ambiguous"

sexcheck["F_LABEL"] = sexcheck["F"].apply(f_label)

merged = sexcheck.merge(sex_info[["IID", "COLLECTED_SEX"]], on="IID", how="left")
merged["COLLECTED_SEX"] = merged["COLLECTED_SEX"].fillna("Unknown")

# Discordant = PLINK called a sex AND it disagrees with collected
def is_discordant(row):
    if row["INFERRED_SEX"] == "Unknown" or row["COLLECTED_SEX"] == "Unknown":
        return False
    return row["INFERRED_SEX"] != row["COLLECTED_SEX"]

merged["DISCORDANT"] = merged.apply(is_discordant, axis=1)

merged.to_csv(args.out_annot, sep="\t", index=False)

discord = merged[merged["DISCORDANT"]][["FID", "IID", "COLLECTED_SEX", "INFERRED_SEX", "F", "STATUS"]]
discord.to_csv(args.out_discord, sep="\t", index=False)

print(f"Total samples checked : {len(merged)}")
print(f"Sex-discordant samples: {merged['DISCORDANT'].sum()}")
print(f"PLINK STATUS=PROBLEM  : {(merged['STATUS'] == 'PROBLEM').sum()}")
print(f"Ambiguous F-statistic : {(merged['F_LABEL'] == 'Ambiguous').sum()}")