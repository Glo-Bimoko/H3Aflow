"""
extract_xy_intensity.py
=======================
Reads the concatenated per-sample GTC TSV produced by GTC_TO_VCF (via
bcftools +gtc2vcf) and extracts per-sample mean X and Y raw intensities.
Merges with collected sex from sex_info and writes:
  --out   xy_intensity.tsv   (per-sample summary table)
  --plot  xy_intensity_plot.png  (scatter: mean_X vs mean_Y, coloured by sex)

Expected TSV columns (bcftools gtc2vcf --format GTC output):
  SAMPLE_ID  CHR  NORMX  NORMY  [other columns ignored]

The script is tolerant of column-name casing and common aliases.
"""

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Extract XY intensities from GTC TSV")
parser.add_argument("--tsv",      required=True,  help="Concatenated GTC TSV (all plates)")
parser.add_argument("--sex_info", required=True,  help="TSV: sampleid, sex (0=Female, 1=Male)")
parser.add_argument("--out",      required=True,  help="Output TSV path")
parser.add_argument("--plot",     required=True,  help="Output PNG path")
args = parser.parse_args()

# ── Load GTC TSV ───────────────────────────────────────────────────────────────
print(f"[extract_xy_intensity] Reading {args.tsv} …", flush=True)
gtc = pd.read_csv(args.tsv, sep="\t", low_memory=False)
gtc.columns = gtc.columns.str.strip().str.upper()

# Resolve column aliases
SAMPLE_ALIASES = ["SAMPLE_ID", "SAMPLEID", "SAMPLE", "ID"]
X_ALIASES      = ["NORMX", "X_NORM", "X_RAW", "NORMALISEDX", "X"]
Y_ALIASES      = ["NORMY", "Y_NORM", "Y_RAW", "NORMALISEDY", "Y"]
CHR_ALIASES    = ["CHR", "CHROM", "#CHROM", "CHROMOSOME"]

def resolve_col(df, aliases, label):
    for a in aliases:
        if a in df.columns:
            return a
    sys.exit(
        f"[extract_xy_intensity] ERROR: Could not find {label} column.\n"
        f"  Tried: {aliases}\n  Found: {list(df.columns)}"
    )

col_sample = resolve_col(gtc, SAMPLE_ALIASES, "SAMPLE_ID")
col_x      = resolve_col(gtc, X_ALIASES,      "NORMX")
col_y      = resolve_col(gtc, Y_ALIASES,       "NORMY")
col_chr    = resolve_col(gtc, CHR_ALIASES,     "CHR")

# Coerce intensity columns to numeric (bcftools may write '.' for missing)
gtc[col_x] = pd.to_numeric(gtc[col_x], errors="coerce")
gtc[col_y] = pd.to_numeric(gtc[col_y], errors="coerce")

# ── Compute per-sample means (all autosomes + sex chrs) ───────────────────────
# Also compute means restricted to chrX and chrY for sex QC context
chrom_str = gtc[col_chr].astype(str).str.upper().str.lstrip("CHR")

autosome_mask = chrom_str.str.match(r"^\d+$")
chrx_mask     = chrom_str.isin(["X", "23"])
chry_mask     = chrom_str.isin(["Y", "24"])

def mean_intensity(sub, x_col, y_col, sample_col):
    """Return per-sample mean X and Y from a subset of rows."""
    return (
        sub.groupby(sample_col)[[x_col, y_col]]
        .mean()
        .rename(columns={x_col: "MEAN_X", y_col: "MEAN_Y"})
    )

df_all  = mean_intensity(gtc,                    col_x, col_y, col_sample)
df_auto = mean_intensity(gtc[autosome_mask],     col_x, col_y, col_sample).add_suffix("_AUTO")
df_x    = mean_intensity(gtc[chrx_mask],         col_x, col_y, col_sample).add_suffix("_X")
df_y    = mean_intensity(gtc[chry_mask],         col_x, col_y, col_sample).add_suffix("_Y")

summary = (
    df_all
    .join(df_auto, how="left")
    .join(df_x,    how="left")
    .join(df_y,    how="left")
    .reset_index()
    .rename(columns={col_sample: "IID"})
)

# ── Merge sex_info ─────────────────────────────────────────────────────────────
print(f"[extract_xy_intensity] Reading {args.sex_info} …", flush=True)
sex_info = pd.read_csv(args.sex_info, sep="\t", dtype=str)
sex_info.columns = sex_info.columns.str.strip().str.lower()
sex_info = sex_info.rename(columns={"sampleid": "IID"})
sex_info["IID"] = sex_info["IID"].str.strip()
sex_info["COLLECTED_SEX"] = sex_info["sex"].map(
    {"0": "Female", "1": "Male"}
).fillna("Unknown")

summary = summary.merge(sex_info[["IID", "COLLECTED_SEX"]], on="IID", how="left")
summary["COLLECTED_SEX"] = summary["COLLECTED_SEX"].fillna("Unknown")

# ── Write TSV ──────────────────────────────────────────────────────────────────
summary.to_csv(args.out, sep="\t", index=False)
print(f"[extract_xy_intensity] Written {len(summary)} samples → {args.out}", flush=True)

# ── Plot ───────────────────────────────────────────────────────────────────────
SEX_PALETTE = {
    "Male":    "#2196F3",   # blue
    "Female":  "#E91E63",   # pink
    "Unknown": "#9E9E9E",   # grey
}

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle("X / Y Raw Intensity QC", fontsize=14, fontweight="bold")

# --- Panel 1: mean chrX vs mean chrY intensity (sex separation plot) ----------
ax = axes[0]
for sex, grp in summary.groupby("COLLECTED_SEX"):
    ax.scatter(
        grp["MEAN_X_X"], grp["MEAN_Y_Y"],
        c=SEX_PALETTE.get(sex, "#9E9E9E"),
        alpha=0.6, s=18, label=sex, edgecolors="none"
    )
ax.set_xlabel("Mean chrX Intensity (NORMX)", fontsize=11)
ax.set_ylabel("Mean chrY Intensity (NORMY)", fontsize=11)
ax.set_title("chrX vs chrY Intensity\n(coloured by collected sex)", fontsize=11)
ax.legend(title="Collected Sex", framealpha=0.7)
ax.grid(True, linewidth=0.4, alpha=0.5)

# --- Panel 2: genome-wide mean X vs mean Y intensity (overall signal level) ---
ax = axes[1]
for sex, grp in summary.groupby("COLLECTED_SEX"):
    ax.scatter(
        grp["MEAN_X"], grp["MEAN_Y"],
        c=SEX_PALETTE.get(sex, "#9E9E9E"),
        alpha=0.6, s=18, label=sex, edgecolors="none"
    )
ax.set_xlabel("Genome-wide Mean X Intensity (NORMX)", fontsize=11)
ax.set_ylabel("Genome-wide Mean Y Intensity (NORMY)", fontsize=11)
ax.set_title("Genome-wide X vs Y Intensity\n(coloured by collected sex)", fontsize=11)
ax.legend(title="Collected Sex", framealpha=0.7)
ax.grid(True, linewidth=0.4, alpha=0.5)

plt.tight_layout()
plt.savefig(args.plot, dpi=150, bbox_inches="tight")
plt.close()
print(f"[extract_xy_intensity] Plot saved → {args.plot}", flush=True)

# ── Summary stats ──────────────────────────────────────────────────────────────
for sex, grp in summary.groupby("COLLECTED_SEX"):
    print(
        f"  {sex:8s}: n={len(grp):5d}  "
        f"chrX mean_X={grp['MEAN_X_X'].mean():.4f}  "
        f"chrY mean_Y={grp['MEAN_Y_Y'].mean():.4f}"
    )
