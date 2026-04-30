"""
annotate_sex_check.py
Merges PLINK --check-sex output with collected gender from the sex_info file.
Flags samples where inferred sex disagrees with collected sex.

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
args = parser.parse_args()

sexcheck = pd.read_csv(args.sexcheck, sep=r"\s+")

sex_info = pd.read_csv(args.sex_info, sep="\t", dtype=str)
sex_info.columns = sex_info.columns.str.strip().str.lower()
sex_info = sex_info.rename(columns={"sampleid": "IID"})
sex_info["IID"] = sex_info["IID"].astype(str).str.strip()

# Map collected sex: 0→Female, 1→Male
sex_info["COLLECTED_SEX"] = sex_info["sex"].map({"0": "Female", "1": "Male"}).fillna("Unknown")

# Map PLINK inferred sex: 1→Male, 2→Female, 0→Unknown
sexcheck["INFERRED_SEX"] = sexcheck["SNPSEX"].map({1: "Male", 2: "Female", 0: "Unknown"})

# Inferred F-statistic interpretation label
def f_label(f):
    if pd.isna(f):    return "Unknown"
    if f > 0.8:       return "Male"
    if f < 0.2:       return "Female"
    return "Ambiguous"
sexcheck["F_LABEL"] = sexcheck["F"].apply(f_label)

merged = sexcheck.merge(sex_info[["IID","COLLECTED_SEX"]], on="IID", how="left")
merged["COLLECTED_SEX"] = merged["COLLECTED_SEX"].fillna("Unknown")

# Discordant = PLINK called a sex AND it disagrees with collected
def is_discordant(row):
    if row["INFERRED_SEX"] == "Unknown" or row["COLLECTED_SEX"] == "Unknown":
        return False
    return row["INFERRED_SEX"] != row["COLLECTED_SEX"]

merged["DISCORDANT"] = merged.apply(is_discordant, axis=1)

merged.to_csv(args.out_annot, sep="\t", index=False)

discord = merged[merged["DISCORDANT"]][["FID","IID","COLLECTED_SEX","INFERRED_SEX","F","STATUS"]]
discord.to_csv(args.out_discord, sep="\t", index=False)

print(f"Total samples checked : {len(merged)}")
print(f"Sex-discordant samples: {merged['DISCORDANT'].sum()}")
print(f"PLINK STATUS=PROBLEM  : {(merged['STATUS']=='PROBLEM').sum()}")
print(f"Ambiguous F-statistic : {(merged['F_LABEL']=='Ambiguous').sum()}")
