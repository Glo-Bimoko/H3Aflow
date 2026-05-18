#!/usr/bin/env python3
"""
extract_xy_intensity.py
=======================
Reads the concatenated per-sample GTC TSV produced by GTC_TO_VCF (via
bcftools +gtc2vcf) and extracts per-sample mean X and Y raw intensities.
Merges with collected sex from sex_info and writes:
  --out   xy_intensity.tsv   (per-sample summary table)
  --plot  xy_intensity_plot.html  (interactive scatter plot using Plotly)

Expected TSV columns (bcftools gtc2vcf --format GTC output):
  SAMPLE_ID  CHR  POS  REF  ALT  NORMX  NORMY

The script is tolerant of column-name casing and common aliases.
Memory-efficient: reads the TSV in chunks so large files (10+ GB) are handled
without loading everything into RAM at once.
"""

import argparse
import sys
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from collections import defaultdict

# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Extract XY intensities from GTC TSV")
parser.add_argument("--tsv",       required=True,  help="Concatenated GTC TSV (all plates)")
parser.add_argument("--sex_info",  required=True,  help="TSV: sampleid, sex (0=Female, 1=Male)")
parser.add_argument("--out",       required=True,  help="Output TSV path")
parser.add_argument("--plot",      required=True,  help="Output HTML path")
parser.add_argument("--chunksize", type=int, default=2_000_000,
                    help="Rows per chunk when reading TSV (default: 2000000)")
args = parser.parse_args()

# ── Column names (for headerless files) ─────────────────────────────────────────
COL_NAMES = ["SAMPLE_ID", "CHR", "POS", "REF", "ALT", "NORMX", "NORMY"]

# ── Function to read TSV (handles headerless files) ────────────────────────────
def read_tsv_smart(filepath, chunksize=None):
    """Read TSV file, detecting if it has headers or not"""
    with open(filepath, 'r') as f:
        first_line = f.readline().strip().split('\t')
    
    # Check if first line contains expected column names
    has_header = any(col.upper() in [c.upper() for c in COL_NAMES] for col in first_line)
    
    if chunksize:
        if has_header:
            return pd.read_csv(filepath, sep="\t", chunksize=chunksize, low_memory=False)
        else:
            return pd.read_csv(filepath, sep="\t", header=None, names=COL_NAMES, 
                             chunksize=chunksize, low_memory=False)
    else:
        if has_header:
            return pd.read_csv(filepath, sep="\t", low_memory=False)
        else:
            return pd.read_csv(filepath, sep="\t", header=None, names=COL_NAMES, low_memory=False)

# ── Read file and detect columns ────────────────────────────────────────────────
print(f"[extract_xy_intensity] Reading {args.tsv} …", flush=True)

# Peek at first chunk to resolve column names
first_chunk = None
for chunk in read_tsv_smart(args.tsv, chunksize=1000):
    first_chunk = chunk
    break

if first_chunk is None:
    sys.exit(f"[extract_xy_intensity] ERROR: Could not read any data from {args.tsv}")

# Ensure column names are uppercase
first_chunk.columns = first_chunk.columns.str.strip().str.upper()
all_cols = list(first_chunk.columns)

# Find required columns (case-insensitive)
col_sample = None
col_chr = None
col_x = None
col_y = None

for col in all_cols:
    col_upper = col.upper()
    if col_sample is None and col_upper in ["SAMPLE_ID", "SAMPLEID", "SAMPLE", "ID"]:
        col_sample = col
    if col_chr is None and col_upper in ["CHR", "CHROM", "#CHROM", "CHROMOSOME"]:
        col_chr = col
    if col_x is None and col_upper in ["NORMX", "X_NORM", "X_RAW", "NORMALISEDX", "X"]:
        col_x = col
    if col_y is None and col_upper in ["NORMY", "Y_NORM", "Y_RAW", "NORMALISEDY", "Y"]:
        col_y = col

if not all([col_sample, col_chr, col_x, col_y]):
    sys.exit(f"[extract_xy_intensity] ERROR: Could not find required columns.\n"
             f"  Tried aliases for SAMPLE_ID, CHR, NORMX, NORMY\n"
             f"  Found: {all_cols}")

needed_cols = [col_sample, col_chr, col_x, col_y]
print(f"[extract_xy_intensity] Using columns: {needed_cols}", flush=True)

# ── Accumulators for chunked mean computation ──────────────────────────────────
# accum[group][sample] = [sum_x, sum_y, count]
accum = {g: defaultdict(lambda: [0.0, 0.0, 0]) for g in ("all", "auto", "X", "Y")}

chunk_n = 0
for chunk in read_tsv_smart(args.tsv, chunksize=args.chunksize):
    # Ensure column names are uppercase
    chunk.columns = chunk.columns.str.strip().str.upper()
    
    chunk[col_x] = pd.to_numeric(chunk[col_x], errors="coerce")
    chunk[col_y] = pd.to_numeric(chunk[col_y], errors="coerce")
    chunk = chunk.dropna(subset=[col_x, col_y])

    # Process chromosome information
    chrom_str = chunk[col_chr].astype(str).str.upper().str.replace('^CHR', '', regex=True)
    auto_mask = chrom_str.str.match(r'^\d+$') & ~chrom_str.str.match(r'^23$|^24$')
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

# ── Create interactive Plotly HTML report ──────────────────────────────────────
print(f"[extract_xy_intensity] Creating interactive plot…", flush=True)

# Color mapping
SEX_COLORS = {
    "Male": "#2196F3",
    "Female": "#E91E63",
    "Unknown": "#9E9E9E"
}

# Create subplots with 2 columns
fig = make_subplots(
    rows=1, cols=2,
    subplot_titles=("chrX vs chrY Intensity", "Genome-wide X vs Y Intensity"),
    horizontal_spacing=0.15
)

# Add traces for first subplot (chrX vs chrY)
for sex, color in SEX_COLORS.items():
    df_sex = summary[summary["COLLECTED_SEX"] == sex]
    if not df_sex.empty:
        fig.add_trace(
            go.Scatter(
                x=df_sex["MEAN_X_X"],
                y=df_sex["MEAN_Y_Y"],
                mode='markers',
                name=sex,
                marker=dict(
                    size=8,
                    color=color,
                    opacity=0.7,
                    line=dict(width=0)
                ),
                text=df_sex["IID"],
                hovertemplate="<b>%{text}</b><br>" +
                              "Mean chrX Intensity: %{x:.3f}<br>" +
                              "Mean chrY Intensity: %{y:.3f}<br>" +
                              "<extra></extra>"
            ),
            row=1, col=1
        )

# Add traces for second subplot (genome-wide)
for sex, color in SEX_COLORS.items():
    df_sex = summary[summary["COLLECTED_SEX"] == sex]
    if not df_sex.empty:
        fig.add_trace(
            go.Scatter(
                x=df_sex["MEAN_X"],
                y=df_sex["MEAN_Y"],
                mode='markers',
                name=sex,
                marker=dict(
                    size=8,
                    color=color,
                    opacity=0.7,
                    line=dict(width=0)
                ),
                text=df_sex["IID"],
                hovertemplate="<b>%{text}</b><br>" +
                              "Genome-wide Mean X: %{x:.3f}<br>" +
                              "Genome-wide Mean Y: %{y:.3f}<br>" +
                              "<extra></extra>",
                showlegend=False  # Hide duplicate legend entries
            ),
            row=1, col=2
        )

# Update layout
fig.update_layout(
    title=dict(
        text="X / Y Raw Intensity QC - Interactive Report",
        font=dict(size=16, weight='bold'),
        x=0.5
    ),
    hovermode='closest',
    width=1200,
    height=600,
    legend=dict(
        title="Collected Sex",
        x=1.02,
        y=1,
        xanchor='left',
        bgcolor='rgba(255, 255, 255, 0.8)',
        bordercolor='black',
        borderwidth=1
    ),
    template='plotly_white'
)

# Update axes
fig.update_xaxes(title_text="Mean chrX Intensity (NORMX)", row=1, col=1, gridcolor='lightgray', gridwidth=0.5)
fig.update_yaxes(title_text="Mean chrY Intensity (NORMY)", row=1, col=1, gridcolor='lightgray', gridwidth=0.5)
fig.update_xaxes(title_text="Genome-wide Mean X Intensity (NORMX)", row=1, col=2, gridcolor='lightgray', gridwidth=0.5)
fig.update_yaxes(title_text="Genome-wide Mean Y Intensity (NORMY)", row=1, col=2, gridcolor='lightgray', gridwidth=0.5)

# Save as HTML
fig.write_html(args.plot)
print(f"[extract_xy_intensity] Interactive plot saved → {args.plot}", flush=True)

# ── Summary stats ──────────────────────────────────────────────────────────────
print("\n[extract_xy_intensity] Summary statistics by sex:")
for sex, grp in summary.groupby("COLLECTED_SEX"):
    print(
        f"  {sex:8s}: n={len(grp):5d}  "
        f"chrX mean_X={grp['MEAN_X_X'].mean():.4f}  "
        f"chrY mean_Y={grp['MEAN_Y_Y'].mean():.4f}"
    )