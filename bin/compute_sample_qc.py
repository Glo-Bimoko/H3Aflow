"""
compute_sample_qc.py
Reads PLINK .imiss and .het files, flags samples failing call rate or
heterozygosity thresholds, and writes pass/fail lists + stats table.
"""
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--imiss",    required=True)
parser.add_argument("--het",      required=True)
parser.add_argument("--mind",     type=float, default=0.05)
parser.add_argument("--het_sd",   type=float, default=3)
parser.add_argument("--out_stats",required=True)
parser.add_argument("--out_pass", required=True)
parser.add_argument("--out_fail", required=True)
args = parser.parse_args()

imiss = pd.read_csv(args.imiss, sep=r"\s+")
het   = pd.read_csv(args.het,   sep=r"\s+")

# Observed heterozygosity: (N(NM) - O(HOM)) / N(NM)
het["OBS_HET"] = (het["N(NM)"] - het["O(HOM)"]) / het["N(NM)"]

df = imiss.merge(het[["FID","IID","OBS_HET"]], on=["FID","IID"])

# Flags
df["FAIL_MIND"] = df["F_MISS"] > args.mind
mean_het = df["OBS_HET"].mean()
sd_het   = df["OBS_HET"].std()
df["FAIL_HET"]  = (df["OBS_HET"] < mean_het - args.het_sd * sd_het) | \
                  (df["OBS_HET"] > mean_het + args.het_sd * sd_het)
df["FAIL_ANY"]  = df["FAIL_MIND"] | df["FAIL_HET"]
df["FAIL_REASON"] = df.apply(
    lambda r: ",".join(filter(None, [
        "HIGH_MISSINGNESS" if r["FAIL_MIND"] else "",
        "HET_OUTLIER"      if r["FAIL_HET"]  else "",
    ])) or "PASS",
    axis=1
)

df.to_csv(args.out_stats, sep="\t", index=False)

pass_df = df[~df["FAIL_ANY"]][["FID","IID"]]
fail_df = df[ df["FAIL_ANY"]][["FID","IID"]]

pass_df.to_csv(args.out_pass, sep="\t", index=False, header=False)
fail_df.to_csv(args.out_fail, sep="\t", index=False, header=False)

print(f"Samples passing QC : {len(pass_df)}")
print(f"Samples failing QC : {len(fail_df)}")
print(f"  High missingness : {df['FAIL_MIND'].sum()}")
print(f"  Het outlier      : {df['FAIL_HET'].sum()}")
