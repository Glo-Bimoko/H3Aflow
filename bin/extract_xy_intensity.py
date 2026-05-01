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
Memory-efficient: reads the TSV in chunks so large files (10+ GB) are handled
without loading everything into RAM at once.
"""

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Extract XY intensities from GTC TSV")
parser.add_argument("--tsv",       required=True,  help="Concatenated GTC TSV (all plates)")
parser.add_argument("--sex_info",  required=True,  help="TSV: sampleid, sex (0=Female, 1=Male)")
parser.add_argument("--out",       required=True,  help="Output TSV path")
parser.add_argument("--plot",      required=True,  help="Output PNG path")
parser.add_argument("--chunksize", type=int, default=2_000_000,
                    help="Rows per chunk when reading TSV (default: 2000000)")
args = parser.parse_args()

# ── Column aliases ─────────────────────────────────────────────────────────────
SAMPLE_ALIASES = ["SAMPLE_ID", "SAMPLEID", "SAMPLE", "ID"]
X_ALIASES      = ["NORMX", "X_NORM", "X_RAW", "NORMALISEDX", "X"]
Y_ALIASES      = ["NORMY", "Y_NORM", "Y_RAW", "NORMALISEDY", "Y"]
CHR_ALIASES    = ["CHR", "CHROM", "#CHROM", "CHROMOSOME"]

def resolve_col(columns, aliases, label):
    for a in aliases:
        if a in columns:
            return a
    sys.exit(
        f"[extract_xy_intensity] ERROR: Could not find {label} column.\n"
        f"  Tried: {aliases}\n  Found: {list(columns)}"
    )

# ── Peek at header to resolve column names and select only needed cols ─────────
print(f"[extract_xy_intensity] Reading {args.tsv} …", flush=True)

header_df = pd.read_csv(args.tsv, sep="\t", nrows=0)
header_df.columns = header_df.columns.str.strip().str.upper()
all_cols = list(header_df.columns)

col_sample = resolve_col(all_cols, SAMPLE_ALIASES, "SAMPLE_ID")
col_x      = resolve_col(all_cols, X_ALIASES,      "NORMX")
col_y      = resolve_col(all_cols, Y_ALIASES,       "NORMY")
col_chr    = resolve_col(all_cols, CHR_ALIASES,     "CHR")

needed_cols = [col_sample, col_chr, col_x, col_y]
print(f"[extract_xy_intensity] Using columns: {needed_cols}", flush=True)

# ── Accumulators for chunked mean computation ──────────────────────────────────
# For each (sample, group) we track sum_x, sum_y, count so we can compute means
# without holding all rows in memory.
# Groups: 'all', 'auto', 'X', 'Y'

from collections import defaultdict

# accum[group][sample] = [sum_x, sum_y, count]
accum = {g: defaultdict(lambda: [0.0, 0.0, 0]) for g in ("all", "auto", "X", "Y")}

chunk_n = 0
for chunk in pd.read_csv(
    args.tsv, sep="\t",
    usecols=needed_cols,
    chunksize=args.chunksize,
    low_memory=False
):
    chunk.columns = chunk.columns.str.strip().str.upper()
    chunk[col_x] = pd.to_numeric(chunk[col_x], errors="coerce")
    chunk[col_y] = pd.to_numeric(chunk[col_y], errors="coerce")
    chunk = chunk.dropna(subset=[col_x, col_y])

    chrom_str = chunk[col_chr].astype(str).str.upper().str.lstrip("CHR")
    auto_mask = chrom_str.str.match(r"^\d+$")
    chrx_mask = chrom_str.isin(["X", "23"])
    chry_mask = chrom_str.isin(["Y", "24"])

    masks = {
        "all":  pd.Series(True, index=chunk.index),
        "auto": auto_mask,
        "X":    chrx_mask,
        "Y":    chry_mask,
    }

    for grp_name, mask in masks.items():
        sub = chunk[mask]
        if sub.empty:
            continue
        for sample, grp in sub.groupby(col_sample):
            a = accum[grp_name][sample]
            a[0] += grp[col_x].sum()
            a[1] += grp[col_y].sum()
            a[2] += len(grp)

    chunk_n += 1
    if chunk_n % 10 == 0:
        print(f"[extract_xy_intensity] Processed {chunk_n * args.chunksize:,} rows …",
              flush=True)

print(f"[extract_xy_intensity] Finished reading. Building summary …", flush=True)

# ── Build summary dataframe from accumulators ──────────────────────────────────
def accum_to_df(acc, x_col, y_col):
    rows = []
    for sample, (sx, sy, n) in acc.items():
        rows.append({
            "IID":  sample,
            x_col: sx / n if n > 0 else np.nan,
            y_col: sy / n if n > 0 else np.nan,
        })
    return pd.DataFrame(rows).set_index("IID") if rows else pd.DataFrame()

df_all  = accum_to_df(accum["all"],  "MEAN_X",      "MEAN_Y")
df_auto = accum_to_df(accum["auto"], "MEAN_X_AUTO", "MEAN_Y_AUTO")
df_x    = accum_to_df(accum["X"],    "MEAN_X_X",    "MEAN_Y_X")
df_y    = accum_to_df(accum["Y"],    "MEAN_X_Y",    "MEAN_Y_Y")

summary = (
    df_all
    .join(df_auto, how="left")
    .join(df_x,    how="left")
    .join(df_y,    how="left")
    .reset_index()
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
    "Male":    "#2196F3",
    "Female":  "#E91E63",
    "Unknown": "#9E9E9E",
}

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle("X / Y Raw Intensity QC", fontsize=14, fontweight="bold")

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