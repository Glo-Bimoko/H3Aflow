#!/usr/bin/env python3
"""
generate_report.py
==================
Compiles a self-contained HTML QC report for the H3aFlow pipeline.
Reads all upstream QC outputs and writes:
  --out_html     qc_report.html        (standalone HTML; no external deps)
  --out_flagged  flagged_samples.tsv   (union of all QC failures)
  --out_summary  cohort_summary.tsv    (one-row key-value table)

Inputs (all from prior pipeline stages):
  --sexcheck      sexcheck.txt               (GTC computed_gender vs collected; gtc_sex_check.py)
  --xy_tsv        xy_intensity.tsv           (median X/Y from GTC stats; gtc_sex_check.py)
  --qc_stats      sample_qc_stats.tsv        (compute_sample_qc.py)
  --genome        ibd.genome                 (PLINK --genome)
  --eigenvec      pca.eigenvec               (PLINK2 --pca)
  --sex_info      sex_info.tsv               (raw collected sex)
  --concordance   pairwise_concordance.tsv   (pairwise_concordance.py) [optional]
  --samplesheet   samplesheet.csv            (original samplesheet with Plate/Well) [optional]
"""

import argparse
import base64
import io
import os
import sys
import textwrap
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Generate HTML QC report")
parser.add_argument("--sexcheck",     required=True)
parser.add_argument("--xy_tsv",      required=True)
parser.add_argument("--qc_stats",    required=True)
parser.add_argument("--genome",      required=True)
parser.add_argument("--eigenvec",    required=True)
parser.add_argument("--sex_info",    required=True)
parser.add_argument("--out_html",    required=True)
parser.add_argument("--out_flagged", required=True)
parser.add_argument("--out_summary", required=True)
# NEW: concordance and samplesheet (optional)
parser.add_argument("--concordance", default=None,
                    help="pairwise_concordance.tsv from pairwise_concordance.py")
parser.add_argument("--samplesheet", default=None,
                    help="Original samplesheet CSV with Plate Number and Well Position columns")
parser.add_argument("--gtc_qc_summary", default=None,
                    help="GTC-level per-sample QC summary TSV (from GTC_QC)")
parser.add_argument("--poorgc10", default=None,
                    help="List of samples failing p10_gc threshold (poorgc10.lst)")
parser.add_argument("--gc10_threshold", type=float, default=0.15,
                    help="Minimum p10_gc / gencall_score_10_percentile (matches nextflow.config)")
# Optional thresholds (match nextflow.config defaults)
parser.add_argument("--pi_hat",              type=float, default=0.1875)
parser.add_argument("--mind",                type=float, default=0.05)
parser.add_argument("--het_sd",              type=float, default=3.0)
parser.add_argument("--concordance_warn",    type=float, default=84.0,
                    help="Concordance %% above which a pair is included in the report (default: 84)")
parser.add_argument("--concordance_flag",    type=float, default=99.0,
                    help="Concordance %% above which a pair is flagged as likely duplicate (default: 99)")
args = parser.parse_args()

# ══════════════════════════════════════════════════════════════════════════════
# 1.  LOAD ALL DATA
# ══════════════════════════════════════════════════════════════════════════════

print("[generate_report] Loading input files …", flush=True)

# ── sex_info ──────────────────────────────────────────────────────────────────
sex_info = pd.read_csv(args.sex_info, sep="\t", dtype=str)
sex_info.columns = sex_info.columns.str.strip().str.lower()
sex_info = sex_info.rename(columns={"sampleid": "IID"})
sex_info["IID"] = sex_info["IID"].str.strip()
sex_info["COLLECTED_SEX"] = sex_info["sex"].map({"0": "Female", "1": "Male"}).fillna("Unknown")

total_samples = len(sex_info)

# ── sample QC stats ───────────────────────────────────────────────────────────
qc = pd.read_csv(args.qc_stats, sep="\t")
qc["IID"] = qc["IID"].astype(str).str.strip()
n_pass_mind  = (~qc["FAIL_MIND"]).sum()
n_fail_mind  = qc["FAIL_MIND"].sum()
n_fail_het   = qc["FAIL_HET"].sum()
n_fail_any   = qc["FAIL_ANY"].sum()
n_pass_any   = (~qc["FAIL_ANY"]).sum()

# ── sexcheck ──────────────────────────────────────────────────────────────────
sexcheck = pd.read_csv(args.sexcheck, sep=r"\s+")
sexcheck["IID"] = sexcheck["IID"].astype(str).str.strip()

if "COLLECTED_SEX" not in sexcheck.columns:
    sexcheck["INFERRED_SEX"] = sexcheck["SNPSEX"].map({1: "Male", 2: "Female", 0: "Unknown"})
    sexcheck = sexcheck.merge(sex_info[["IID","COLLECTED_SEX"]], on="IID", how="left")
    sexcheck["COLLECTED_SEX"] = sexcheck["COLLECTED_SEX"].fillna("Unknown")
    def _is_discord(row):
        if row["INFERRED_SEX"] == "Unknown" or row["COLLECTED_SEX"] == "Unknown":
            return False
        return row["INFERRED_SEX"] != row["COLLECTED_SEX"]
    sexcheck["DISCORDANT"] = sexcheck.apply(_is_discord, axis=1)

n_sex_discord = int(sexcheck["DISCORDANT"].sum()) if "DISCORDANT" in sexcheck.columns else "N/A"
n_sex_problem = int((sexcheck["STATUS"] == "PROBLEM").sum()) if "STATUS" in sexcheck.columns else "N/A"

# ── XY intensities ────────────────────────────────────────────────────────────
xy = pd.read_csv(args.xy_tsv, sep="\t")
xy["IID"] = xy["IID"].astype(str).str.strip()

# ── IBD genome ────────────────────────────────────────────────────────────────
genome = pd.read_csv(args.genome, sep=r"\s+")
genome.columns = genome.columns.str.strip()
flagged_ibd = genome[genome["PI_HAT"] >= args.pi_hat] if "PI_HAT" in genome.columns else pd.DataFrame()

def _rel(pi):
    if pi >= 0.90:  return "Duplicate/MZ_twin"
    if pi >= 0.45:  return "1st_degree"
    if pi >= 0.20:  return "2nd_degree"
    return "Cryptic"

if len(flagged_ibd):
    flagged_ibd = flagged_ibd.copy()
    flagged_ibd["RELATIONSHIP"] = flagged_ibd["PI_HAT"].apply(_rel)

n_ibd_pairs    = len(flagged_ibd)
n_ibd_dup      = int((flagged_ibd["RELATIONSHIP"] == "Duplicate/MZ_twin").sum()) if n_ibd_pairs else 0

# ── PCA eigenvec ──────────────────────────────────────────────────────────────
evec = pd.read_csv(args.eigenvec, sep=r"\s+", comment=None)
evec.columns = evec.columns.str.lstrip("#")
evec["IID"] = evec["IID"].astype(str).str.strip()

for col in ["PC1","PC2"]:
    if col in evec.columns:
        m, s = evec[col].mean(), evec[col].std()
        evec[f"{col}_out"] = (evec[col].abs() > m + 6*s) | (evec[col] < m - 6*s)
evec["PCA_OUTLIER"] = (
    evec.get("PC1_out", pd.Series(False, index=evec.index)) |
    evec.get("PC2_out", pd.Series(False, index=evec.index))
)
n_pca_outliers = int(evec["PCA_OUTLIER"].sum())

# ── Samplesheet (optional) ────────────────────────────────────────────────────
samplesheet = None
has_plate_info = False
if args.samplesheet and os.path.exists(args.samplesheet):
    try:
        samplesheet = pd.read_csv(args.samplesheet)
        samplesheet.columns = samplesheet.columns.str.strip().str.lower().str.replace(" ", "_")
        # Normalise common column name variants
        col_renames = {}
        for c in samplesheet.columns:
            if "sample" in c and "id" in c:       col_renames[c] = "sample_id"
            elif "plate" in c:                     col_renames[c] = "plate"
            elif "well" in c:                      col_renames[c] = "well"
            elif "barcode" in c:                   col_renames[c] = "barcode"
        samplesheet = samplesheet.rename(columns=col_renames)
        samplesheet["sample_id"] = samplesheet["sample_id"].astype(str).str.strip()
        has_plate_info = ("plate" in samplesheet.columns and "well" in samplesheet.columns)
        print(f"[generate_report] Samplesheet loaded: {len(samplesheet)} rows, plate_info={has_plate_info}", flush=True)
    except Exception as e:
        print(f"[generate_report] Warning: could not load samplesheet: {e}", flush=True)

# ── Concordance (optional) ────────────────────────────────────────────────────
concordance_df = None
concordance_warn_pairs = None
concordance_flag_pairs = None
n_warn_pairs = 0
n_flag_pairs = 0

if args.concordance and os.path.exists(args.concordance):
    try:
        concordance_df = pd.read_csv(args.concordance, sep="\t")
        concordance_df["SAMPLE_A"] = concordance_df["SAMPLE_A"].astype(str).str.strip()
        concordance_df["SAMPLE_B"] = concordance_df["SAMPLE_B"].astype(str).str.strip()
        # Drop pairs where concordance is NaN
        concordance_df = concordance_df.dropna(subset=["CONCORDANCE_PCT"])
        # Pairs above the warning threshold (84%)
        concordance_warn_pairs = concordance_df[
            concordance_df["CONCORDANCE_PCT"] >= args.concordance_warn
        ].sort_values("CONCORDANCE_PCT", ascending=False).copy()
        # Pairs above the duplicate-flag threshold (99%)
        concordance_flag_pairs = concordance_df[
            concordance_df["CONCORDANCE_PCT"] >= args.concordance_flag
        ].sort_values("CONCORDANCE_PCT", ascending=False).copy()
        n_warn_pairs = len(concordance_warn_pairs)
        n_flag_pairs = len(concordance_flag_pairs)
        print(f"[generate_report] Concordance: {len(concordance_df):,} pairs loaded, "
              f"{n_warn_pairs} ≥{args.concordance_warn}%, {n_flag_pairs} ≥{args.concordance_flag}%", flush=True)
    except Exception as e:
        print(f"[generate_report] Warning: could not load concordance file: {e}", flush=True)

# ── GTC-level QC (optional)
gtc_qc = None
n_poor_gc10 = 0
n_poor_cr = 0
poor_gc10_samples = []
if args.gtc_qc_summary and os.path.exists(args.gtc_qc_summary):
    try:
        gtc_qc = pd.read_csv(args.gtc_qc_summary, sep="\t")
        # Normalise sample id column
        if 'sample_id' in gtc_qc.columns:
            gtc_qc['IID'] = gtc_qc['sample_id'].astype(str).str.strip()
        elif 'IID' in gtc_qc.columns:
            gtc_qc['IID'] = gtc_qc['IID'].astype(str).str.strip()

        if 'pass_gc10' in gtc_qc.columns:
            n_poor_gc10 = int((~gtc_qc['pass_gc10']).sum())
        if 'pass_cr' in gtc_qc.columns:
            n_poor_cr = int((~gtc_qc['pass_cr']).sum())
        print(f"[generate_report] GTC QC summary loaded: {len(gtc_qc)} samples; poor_gc10={n_poor_gc10}, poor_cr={n_poor_cr}", flush=True)
    except Exception as e:
        print(f"[generate_report] Warning: could not load GTC QC summary: {e}", flush=True)

if args.poorgc10 and os.path.exists(args.poorgc10):
    try:
        with open(args.poorgc10) as fh:
            for line in fh:
                parts = line.strip().split()
                if parts:
                    poor_gc10_samples.append(parts[0])
        # If count not determined from summary, set from list
        if n_poor_gc10 == 0:
            n_poor_gc10 = len(poor_gc10_samples)
        print(f"[generate_report] poorgc10 list loaded: {len(poor_gc10_samples)} samples", flush=True)
    except Exception as e:
        print(f"[generate_report] Warning: could not load poorgc10 list: {e}", flush=True)

# ══════════════════════════════════════════════════════════════════════════════
# 2.  BUILD FLAGGED SAMPLES TABLE
# ══════════════════════════════════════════════════════════════════════════════

print("[generate_report] Building flagged samples list …", flush=True)

flag_records = []

for _, row in qc[qc["FAIL_ANY"]].iterrows():
    flag_records.append({"IID": str(row["IID"]), "FLAG": row["FAIL_REASON"], "SOURCE": "SAMPLE_QC"})

if "DISCORDANT" in sexcheck.columns:
    for _, row in sexcheck[sexcheck["DISCORDANT"]].iterrows():
        flag_records.append({"IID": str(row["IID"]), "FLAG": "SEX_DISCORDANT", "SOURCE": "SEX_CHECK"})

if len(flagged_ibd) and "IID2" in flagged_ibd.columns:
    for _, row in flagged_ibd.iterrows():
        flag_records.append({
            "IID": str(row["IID2"]),
            "FLAG": f"IBD_{row['RELATIONSHIP']}",
            "SOURCE": "IBD",
        })

for iid in evec.loc[evec["PCA_OUTLIER"], "IID"].tolist():
    flag_records.append({"IID": str(iid), "FLAG": "PCA_OUTLIER", "SOURCE": "PCA"})

# Concordance duplicates → flag both samples in each pair
if concordance_flag_pairs is not None and len(concordance_flag_pairs):
    flagged_conc_iids = set()
    for _, row in concordance_flag_pairs.iterrows():
        flagged_conc_iids.add(str(row["SAMPLE_A"]))
        flagged_conc_iids.add(str(row["SAMPLE_B"]))
    for iid in flagged_conc_iids:
        flag_records.append({"IID": iid, "FLAG": "CONCORDANCE_DUPLICATE", "SOURCE": "CONCORDANCE"})

# Add GTC-level poor-quality flags (poorgc10 / GTC QC)
if len(poor_gc10_samples):
    for iid in poor_gc10_samples:
        flag_records.append({"IID": str(iid), "FLAG": "POOR_GC10", "SOURCE": "GTC_QC"})
elif gtc_qc is not None and 'pass_gc10' in gtc_qc.columns:
    for _, row in gtc_qc[~gtc_qc['pass_gc10']].iterrows():
        flag_records.append({"IID": str(row['IID']), "FLAG": "POOR_GC10", "SOURCE": "GTC_QC"})

flagged_df = pd.DataFrame(flag_records) if flag_records else pd.DataFrame(columns=["IID","FLAG","SOURCE"])

if len(flagged_df):
    flagged_agg = (
        flagged_df.groupby("IID")
        .agg(FLAGS=("FLAG", lambda x: ",".join(sorted(set(x)))),
             SOURCES=("SOURCE", lambda x: ",".join(sorted(set(x)))))
        .reset_index()
    )
else:
    flagged_agg = pd.DataFrame(columns=["IID","FLAGS","SOURCES"])

flagged_agg.to_csv(args.out_flagged, sep="\t", index=False)
print(f"[generate_report] Flagged samples: {len(flagged_agg)} → {args.out_flagged}", flush=True)

# ══════════════════════════════════════════════════════════════════════════════
# 3.  COHORT SUMMARY TABLE
# ══════════════════════════════════════════════════════════════════════════════

summary_rows = [
    ("Total samples",                  total_samples),
    ("Samples passing call-rate QC",   int(n_pass_mind)),
    ("Samples failing call-rate QC",   int(n_fail_mind)),
    ("Samples failing heterozygosity QC", int(n_fail_het)),
    ("Samples failing ANY QC",         int(n_fail_any)),
    ("Sex-discordant samples",         n_sex_discord),
    ("PLINK sex STATUS=PROBLEM",       n_sex_problem),
    ("IBD pairs flagged",              n_ibd_pairs),
    ("  of which duplicates/MZ",       n_ibd_dup),
    ("PCA outliers (>6 SD)",           n_pca_outliers),
    (f"Concordance pairs ≥{args.concordance_warn}%", n_warn_pairs),
    (f"Concordance pairs ≥{args.concordance_flag}% (likely duplicates)", n_flag_pairs),
    ("Total uniquely flagged samples", len(flagged_agg)),
    ("Samples failing GTC p10_gc (poor cluster quality)", int(n_poor_gc10)),
    ("Samples failing GTC call-rate threshold", int(n_poor_cr)),
]

summary_df = pd.DataFrame(summary_rows, columns=["Metric", "Value"])
summary_df.to_csv(args.out_summary, sep="\t", index=False)
print(f"[generate_report] Summary → {args.out_summary}", flush=True)

# ══════════════════════════════════════════════════════════════════════════════
# 4.  GENERATE EMBEDDED PLOTS
# ══════════════════════════════════════════════════════════════════════════════

def fig_to_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")

SEX_PALETTE = {"Male": "#2196F3", "Female": "#E91E63", "Unknown": "#9E9E9E"}

# ── Plot A: Sample QC ─────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
fig.suptitle("Sample QC – Call Rate & Heterozygosity", fontsize=12, fontweight="bold")

ax = axes[0]
ax.hist(qc["F_MISS"], bins=60, color="#1976D2", edgecolor="white", linewidth=0.4)
ax.axvline(args.mind, color="red", linestyle="--", linewidth=1.3,
           label=f"Cutoff ({args.mind})")
ax.set_xlabel("Missingness rate (F_MISS)", fontsize=10)
ax.set_ylabel("Samples", fontsize=10)
ax.set_title("Call-rate distribution", fontsize=10)
ax.legend(fontsize=8)
ax.grid(True, linewidth=0.4, alpha=0.5)

ax = axes[1]
het_mean = qc["OBS_HET"].mean()
het_sd   = qc["OBS_HET"].std()
ax.hist(qc["OBS_HET"], bins=60, color="#43A047", edgecolor="white", linewidth=0.4)
for sign, ls in [(-1, "--"), (1, "--")]:
    ax.axvline(het_mean + sign * args.het_sd * het_sd,
               color="red", linestyle=ls, linewidth=1.3,
               label=f"±{args.het_sd} SD" if sign == 1 else None)
ax.set_xlabel("Observed heterozygosity", fontsize=10)
ax.set_ylabel("Samples", fontsize=10)
ax.set_title("Heterozygosity distribution", fontsize=10)
ax.legend(fontsize=8)
ax.grid(True, linewidth=0.4, alpha=0.5)

plt.tight_layout()
plot_sampleqc_b64 = fig_to_b64(fig)

# ══════════════════════════════════════════════════════════════════════════
# Per-plate Call Rate vs 10% GenCall scatter plots (GTC-level QC)
# ══════════════════════════════════════════════════════════════════════════
plate_gc_b64s = {}
if gtc_qc is not None:
    merged = gtc_qc.copy()
    # ensure IID and metric columns exist
    if 'IID' in merged.columns:
        merged['IID'] = merged['IID'].astype(str).str.strip()

    # determine call_rate and p10_gc column names with common fallbacks
    call_col = None
    p10_col = None
    for c in ['call_rate','Call_Rate','CALL_RATE','call rate']:
        if c in merged.columns:
            call_col = c
            break
    for c in [
        'p10_gc', 'gencall_score_10_percentile',
        '10%_GC_Score', '10%_GC_SCORE', 'p10gc',
    ]:
        if c in merged.columns:
            p10_col = c
            break

    # attach plate label from samplesheet if available
    if samplesheet is not None and 'sample_id' in samplesheet.columns and 'plate' in samplesheet.columns:
        lookup_plate = samplesheet.set_index('sample_id')['plate'].astype(str).to_dict()
        merged['plate_label'] = merged['IID'].map(lookup_plate).fillna('')
    else:
        # try common plate-like columns in gtc_qc
        plate_col = None
        for c in merged.columns:
            if c.lower().replace('_',' ').replace('%','').strip() in ('institute plate label','institute plate','plate','plate label'):
                plate_col = c
                break
        if plate_col:
            merged['plate_label'] = merged[plate_col].astype(str).str.strip()

    if 'plate_label' in merged.columns and merged['plate_label'].replace('', pd.NA).notna().any():
        for plate_name, grp in merged.groupby('plate_label'):
            try:
                x = grp[p10_col] if p10_col and p10_col in grp.columns else None
                y = grp[call_col] if call_col and call_col in grp.columns else None
                if x is None or y is None:
                    continue
                fig, ax = plt.subplots(figsize=(6,5))
                sexes = grp.get('computed_gender') if 'computed_gender' in grp.columns else None
                colors = None
                if sexes is not None:
                    sex_map = {"0":"#E91E63","1":"#2196F3", "Female":"#E91E63","Male":"#2196F3"}
                    colors = [sex_map.get(str(s), '#9E9E9E') for s in sexes]
                ax.scatter(x, y, c=colors if colors is not None else '#1976D2', alpha=0.6, s=18, edgecolors='none')
                ax.axvline(args.gc10_threshold, color='#C62828', linestyle='--', linewidth=1.2,
                           label=f'p10_gc cutoff ({args.gc10_threshold})')
                ax.axhline(0.95, color='#F57F17', linestyle=':', linewidth=1.0,
                           label='call rate 0.95')
                ax.set_xlabel('10% GenCall score (p10_gc)', fontsize=10)
                ax.set_ylabel('Call rate', fontsize=10)
                ax.set_title(f'Plate: {plate_name} (n={len(grp)})', fontsize=11)
                ax.legend(fontsize=7, loc='lower right')
                ax.grid(True, linewidth=0.4, alpha=0.5)
                ax.set_xlim(left=0)
                ax.set_ylim(0,1.02)
                plate_gc_b64s[str(plate_name)] = fig_to_b64(fig)
            except Exception:
                continue

# ── Plot B: Sex check ─────────────────────────────────────────────────────────
use_gtc_sex = (
    "INFERENCE_METHOD" in sexcheck.columns
    and (sexcheck["INFERENCE_METHOD"] == "GTC_computed_gender").any()
)
sex_plot_title = "Sex Check – GTC computed gender" if use_gtc_sex else "Sex Check – chrX F-Statistic"

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
fig.suptitle(sex_plot_title, fontsize=12, fontweight="bold")

ax = axes[0]
if use_gtc_sex and sexcheck["F"].notna().any():
    for sex, grp in sexcheck.groupby("COLLECTED_SEX"):
        ax.hist(grp["F"].dropna(), bins=50, alpha=0.65,
                color=SEX_PALETTE.get(sex, "#9E9E9E"), label=sex,
                edgecolor="white", linewidth=0.3)
    ax.set_xlabel("log R deviation (GTC)", fontsize=10)
    ax.set_ylabel("Samples", fontsize=10)
    ax.set_title("logR deviation by collected sex", fontsize=10)
elif not use_gtc_sex:
    for sex, grp in sexcheck.groupby("COLLECTED_SEX"):
        ax.hist(grp["F"].dropna(), bins=50, alpha=0.65,
                color=SEX_PALETTE.get(sex, "#9E9E9E"), label=sex,
                edgecolor="white", linewidth=0.3)
    ax.axvline(0.2, color="grey", linestyle=":", linewidth=1)
    ax.axvline(0.8, color="grey", linestyle=":", linewidth=1)
    ax.set_xlabel("F-statistic (chrX inbreeding coefficient)", fontsize=10)
    ax.set_ylabel("Samples", fontsize=10)
    ax.set_title("F-statistic by collected sex", fontsize=10)
else:
    ax.set_visible(False)
if ax.get_visible():
    ax.legend(fontsize=8)
    ax.grid(True, linewidth=0.4, alpha=0.5)

ax = axes[1]
if "MEAN_X_X" in xy.columns and "MEAN_Y_Y" in xy.columns:
    for sex, grp in xy.groupby("COLLECTED_SEX"):
        ax.scatter(grp["MEAN_X_X"], grp["MEAN_Y_Y"],
                   c=SEX_PALETTE.get(sex, "#9E9E9E"),
                   alpha=0.55, s=12, edgecolors="none", label=sex)
    ax.set_xlabel("Mean chrX Intensity", fontsize=10)
    ax.set_ylabel("Mean chrY Intensity", fontsize=10)
    ax.set_title("X/Y Intensity by collected sex", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(True, linewidth=0.4, alpha=0.5)
else:
    ax.set_visible(False)

plt.tight_layout()
plot_sexcheck_b64 = fig_to_b64(fig)

# ── Plot C: IBD ───────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
fig.suptitle("IBD / Relatedness", fontsize=12, fontweight="bold")

ax = axes[0]
if "PI_HAT" in genome.columns and len(genome):
    ax.hist(genome["PI_HAT"], bins=60, color="#7B1FA2",
            edgecolor="white", linewidth=0.3)
    ax.axvline(args.pi_hat, color="red", linestyle="--", linewidth=1.3,
               label=f"Threshold ({args.pi_hat})")
    ax.set_xlabel("PI_HAT", fontsize=10)
    ax.set_ylabel("Pairs", fontsize=10)
    ax.set_title("PI_HAT distribution (all pairs)", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(True, linewidth=0.4, alpha=0.5)
else:
    ax.text(0.5, 0.5, "No pairs in .genome file", ha="center", va="center",
            transform=ax.transAxes, fontsize=11, color="grey")

ax = axes[1]
if len(flagged_ibd) and "Z0" in flagged_ibd.columns and "Z1" in flagged_ibd.columns:
    REL_PAL = {"Duplicate/MZ_twin": "#D32F2F", "1st_degree": "#FF9800",
               "2nd_degree": "#1976D2", "Cryptic": "#616161"}
    for rel, grp in flagged_ibd.groupby("RELATIONSHIP"):
        ax.scatter(grp["Z0"], grp["Z1"], c=REL_PAL.get(rel,"#9E9E9E"),
                   alpha=0.65, s=20, edgecolors="none", label=rel)
    ax.set_xlabel("Z0  (P(IBD=0))", fontsize=10)
    ax.set_ylabel("Z1  (P(IBD=1))", fontsize=10)
    ax.set_title("IBD state proportions (flagged pairs)", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(True, linewidth=0.4, alpha=0.5)
else:
    ax.text(0.5, 0.5, "No flagged IBD pairs", ha="center", va="center",
            transform=ax.transAxes, fontsize=11, color="grey")

plt.tight_layout()
plot_ibd_b64 = fig_to_b64(fig)

# ── Plot D: PCA ───────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle("PCA – Population Structure", fontsize=12, fontweight="bold")

for panel_idx, (pc_a, pc_b) in enumerate([("PC1","PC2"), ("PC1","PC3")]):
    ax = axes[panel_idx]
    if pc_a not in evec.columns or pc_b not in evec.columns:
        ax.set_visible(False)
        continue
    evec_sex = evec.merge(sex_info[["IID","COLLECTED_SEX"]], on="IID", how="left")
    evec_sex["COLLECTED_SEX"] = evec_sex["COLLECTED_SEX"].fillna("Unknown")
    for sex, grp in evec_sex.groupby("COLLECTED_SEX"):
        ax.scatter(grp[pc_a], grp[pc_b],
                   c=SEX_PALETTE.get(sex, "#9E9E9E"),
                   alpha=0.5, s=10, edgecolors="none", label=sex)
    outs = evec_sex[evec_sex["PCA_OUTLIER"]]
    if len(outs):
        ax.scatter(outs[pc_a], outs[pc_b], c="black", marker="x",
                   s=35, linewidths=0.8, label="PCA outlier")
    ax.set_xlabel(pc_a, fontsize=10)
    ax.set_ylabel(pc_b, fontsize=10)
    ax.set_title(f"{pc_a} vs {pc_b}", fontsize=10)
    ax.legend(fontsize=7, framealpha=0.7)
    ax.grid(True, linewidth=0.4, alpha=0.5)

plt.tight_layout()
plot_pca_b64 = fig_to_b64(fig)

# ── Plot E: Concordance distribution ─────────────────────────────────────────
plot_concordance_b64 = None
plate_layout_b64s = {}   # plate_name → b64 png

if concordance_df is not None and len(concordance_df):
    fig, ax = plt.subplots(figsize=(10, 4))
    fig.suptitle("Pairwise Genotype Concordance Distribution", fontsize=12, fontweight="bold")

    ax.hist(concordance_df["CONCORDANCE_PCT"], bins=80, color="#455A64",
            edgecolor="white", linewidth=0.3)
    ax.axvline(args.concordance_warn, color="#F57F17", linestyle="--", linewidth=1.5,
               label=f"Report threshold ({args.concordance_warn}%)")
    ax.axvline(args.concordance_flag, color="#C62828", linestyle="--", linewidth=1.5,
               label=f"Duplicate flag ({args.concordance_flag}%)")
    ax.set_xlabel("Concordance (%)", fontsize=10)
    ax.set_ylabel("Pairs", fontsize=10)
    ax.set_title("All pairwise concordance values", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, linewidth=0.4, alpha=0.5)
    plt.tight_layout()
    plot_concordance_b64 = fig_to_b64(fig)

# ── Plot F: Plate layout heatmaps ────────────────────────────────────────────
def well_to_row_col(well_str):
    """Convert well like 'A01' or 'A1' to (row_idx 0-7, col_idx 0-11) for a 96-well plate."""
    well_str = str(well_str).strip().upper()
    if not well_str or well_str in ("NAN", ""):
        return None, None
    row_letter = well_str[0]
    try:
        col_num = int(well_str[1:])
    except ValueError:
        return None, None
    row_idx = ord(row_letter) - ord('A')
    col_idx = col_num - 1
    if 0 <= row_idx < 8 and 0 <= col_idx < 12:
        return row_idx, col_idx
    return None, None

if has_plate_info and samplesheet is not None:
    # ── Per-sample sex status lookup ──────────────────────────────────────────
    # Keys: sample IID (str); Values: "OK" | "DISCORDANT" | "AMBIGUOUS" | "UNKNOWN"
    sample_sex_status: dict = {}
    if "DISCORDANT" in sexcheck.columns and "STATUS" in sexcheck.columns:
        def _sex_status(row):
            if row.get("DISCORDANT", False):
                return "DISCORDANT"
            inferred = row.get("INFERRED_SEX", "Unknown")
            if inferred == "Unknown":
                return "AMBIGUOUS"
            return "OK"
        for _, srow in sexcheck.iterrows():
            sample_sex_status[str(srow["IID"])] = _sex_status(srow)

    # Build a per-sample max concordance score (highest similarity to any other sample)
    # Only considering pairs above the warn threshold
    sample_max_conc = {}
    sample_conc_partner = {}
    if concordance_warn_pairs is not None and len(concordance_warn_pairs):
        for _, row in concordance_warn_pairs.iterrows():
            sa, sb, pct = str(row["SAMPLE_A"]), str(row["SAMPLE_B"]), row["CONCORDANCE_PCT"]
            if sa not in sample_max_conc or pct > sample_max_conc[sa]:
                sample_max_conc[sa] = pct
                sample_conc_partner[sa] = sb
            if sb not in sample_max_conc or pct > sample_max_conc[sb]:
                sample_max_conc[sb] = pct
                sample_conc_partner[sb] = sa

    plates = samplesheet["plate"].dropna().unique()
    for plate_name in sorted(plates):
        plate_samples = samplesheet[samplesheet["plate"] == plate_name].copy()
        plate_samples["sample_id"] = plate_samples["sample_id"].astype(str).str.strip()

        # 96-well grid: rows A-H (0-7), cols 1-12 (0-11)
        grid_conc  = np.full((8, 12), np.nan)
        grid_label = np.full((8, 12), "", dtype=object)
        grid_is_dup = np.zeros((8, 12), dtype=bool)

        # Sex-discordance status per grid cell
        grid_sex_status = np.full((8, 12), "", dtype=object)

        for _, srow in plate_samples.iterrows():
            r, c = well_to_row_col(srow.get("well", ""))
            if r is None:
                continue
            sid = str(srow["sample_id"])
            max_pct = sample_max_conc.get(sid, np.nan)
            grid_conc[r, c]   = max_pct if not np.isnan(max_pct) else 0.0
            grid_label[r, c]  = sid
            grid_is_dup[r, c] = max_pct >= args.concordance_flag if not np.isnan(max_pct) else False
            grid_sex_status[r, c] = sample_sex_status.get(sid, "UNKNOWN")

        fig, ax = plt.subplots(figsize=(14, 6))
        fig.suptitle(
            f"Plate Layout – {plate_name}\n"
            f"Background: max pairwise concordance per well  |  "
            f"Corner triangle: sex check status",
            fontsize=11, fontweight="bold",
        )

        # Background: grey for empty, colour scale for occupied
        cmap = plt.cm.YlOrRd
        norm = Normalize(vmin=args.concordance_warn, vmax=100)

        # Sex status → triangle fill colour (top-right corner triangle)
        SEX_TRIANGLE_COLOR = {
            "OK":         "#4CAF50",   # green  – concordant
            "DISCORDANT": "#F44336",   # red    – sex discordant
            "AMBIGUOUS":  "#FF9800",   # amber  – ambiguous F-stat
            "UNKNOWN":    "#9E9E9E",   # grey   – no sex-check data
            "":           "#9E9E9E",
        }

        for r in range(8):
            for c in range(12):
                row_letter = chr(ord('A') + r)
                col_num    = c + 1
                well_id    = f"{row_letter}{col_num:02d}"

                # Check if this well is occupied
                match = plate_samples[
                    plate_samples["well"].astype(str).str.strip().str.upper() == well_id
                ]
                if len(match) == 0:
                    # Empty well
                    rect = mpatches.FancyBboxPatch(
                        (c+0.05, r+0.05), 0.9, 0.9,
                        boxstyle="round,pad=0.05",
                        facecolor="#ECEFF1", edgecolor="#CFD8DC", linewidth=0.5,
                    )
                    ax.add_patch(rect)
                else:
                    pct    = grid_conc[r, c]
                    is_dup = grid_is_dup[r, c]
                    sex_st = grid_sex_status[r, c]

                    # ── Background fill (concordance) ─────────────────────────
                    if pct >= args.concordance_warn:
                        face_color = cmap(norm(pct))
                        edge_color = "#B71C1C" if is_dup else "#78909C"
                        edge_lw    = 2.5 if is_dup else 0.8
                    else:
                        face_color = "#E8F5E9"
                        edge_color = "#A5D6A7"
                        edge_lw    = 0.5

                    rect = mpatches.FancyBboxPatch(
                        (c+0.05, r+0.05), 0.9, 0.9,
                        boxstyle="round,pad=0.05",
                        facecolor=face_color, edgecolor=edge_color, linewidth=edge_lw,
                    )
                    ax.add_patch(rect)

                    # ── Sex status triangle (top-right corner) ────────────────
                    # Triangle vertices: top-right corner occupying ~25% of cell
                    tri_size = 0.30
                    tri_x = c + 0.95
                    tri_y = r + 0.05
                    tri = plt.Polygon(
                        [
                            (tri_x,            tri_y),
                            (tri_x - tri_size, tri_y),
                            (tri_x,            tri_y + tri_size),
                        ],
                        closed=True,
                        facecolor=SEX_TRIANGLE_COLOR.get(sex_st, "#9E9E9E"),
                        edgecolor="white",
                        linewidth=0.4,
                        zorder=3,
                    )
                    ax.add_patch(tri)

                    # ── Text labels ───────────────────────────────────────────
                    sid        = grid_label[r, c]
                    label_text = sid[-6:] if len(sid) > 6 else sid
                    pct_text   = f"{pct:.1f}%" if pct >= args.concordance_warn else ""
                    ax.text(c + 0.5, r + 0.62, label_text,
                            ha="center", va="center", fontsize=5.5,
                            color="#212121", fontweight="bold" if is_dup else "normal")
                    if pct_text:
                        ax.text(c + 0.5, r + 0.32, pct_text,
                                ha="center", va="center", fontsize=5,
                                color="#B71C1C" if is_dup else "#5D4037")

        # Axes
        ax.set_xlim(0, 12)
        ax.set_ylim(0, 8)
        ax.set_xticks([c + 0.5 for c in range(12)])
        ax.set_xticklabels([str(i+1) for i in range(12)], fontsize=9)
        ax.set_yticks([r + 0.5 for r in range(8)])
        ax.set_yticklabels([chr(ord('A')+r) for r in range(8)], fontsize=9)
        ax.invert_yaxis()
        ax.set_aspect("equal")
        ax.set_xlabel("Column", fontsize=10)
        ax.set_ylabel("Row", fontsize=10)

        # Colourbar
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.8, pad=0.02)
        cbar.set_label("Max concordance with any other sample (%)", fontsize=9)

        # Legend — two groups: concordance (patches) + sex status (triangles)
        clean_patch  = mpatches.Patch(facecolor="#E8F5E9", edgecolor="#A5D6A7",
                                      label=f"< {args.concordance_warn}% (clean)")
        warn_patch   = mpatches.Patch(facecolor=cmap(norm(90)), edgecolor="#78909C",
                                      label=f"≥ {args.concordance_warn}% (warn)")
        dup_patch    = mpatches.Patch(facecolor=cmap(norm(99.5)), edgecolor="#B71C1C",
                                      linewidth=2.5, label=f"≥ {args.concordance_flag}% (DUPLICATE)")
        empty_patch  = mpatches.Patch(facecolor="#ECEFF1", edgecolor="#CFD8DC",
                                      label="Empty well")
        # Sex triangle legend entries (small squares stand in for triangles)
        sex_ok_p    = mpatches.Patch(facecolor="#4CAF50", edgecolor="white",
                                     label="▲ Sex OK")
        sex_disc_p  = mpatches.Patch(facecolor="#F44336", edgecolor="white",
                                     label="▲ Sex DISCORDANT")
        sex_amb_p   = mpatches.Patch(facecolor="#FF9800", edgecolor="white",
                                     label="▲ Sex AMBIGUOUS")
        sex_unk_p   = mpatches.Patch(facecolor="#9E9E9E", edgecolor="white",
                                     label="▲ Sex UNKNOWN")

        have_sex_data = bool(sample_sex_status)
        sex_legend_handles = (
            [sex_ok_p, sex_disc_p, sex_amb_p, sex_unk_p] if have_sex_data else []
        )
        ax.legend(
            handles=[clean_patch, warn_patch, dup_patch, empty_patch] + sex_legend_handles,
            loc="upper right", fontsize=7.5, bbox_to_anchor=(1.42, 1.0),
            framealpha=0.9,
        )

        plt.tight_layout()
        plate_layout_b64s[plate_name] = fig_to_b64(fig)

# ══════════════════════════════════════════════════════════════════════════════
# 5.  BUILD HTML HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

print("[generate_report] Building HTML …", flush=True)

def summary_table_html(df):
    rows_html = ""
    for _, row in df.iterrows():
        val = row["Value"]
        style = ""
        if isinstance(val, (int, float)) and val > 0 and "passing" not in str(row["Metric"]).lower():
            style = ' style="color:#C62828; font-weight:600;"'
        rows_html += f"<tr><td>{row['Metric']}</td><td{style}>{val}</td></tr>\n"
    return f"""
    <table class="summary-table">
      <thead><tr><th>Metric</th><th>Value</th></tr></thead>
      <tbody>{rows_html}</tbody>
    </table>"""

def flagged_table_html(df, max_rows=500):
    if len(df) == 0:
        return "<p class='ok'>✓ No samples flagged across all QC stages.</p>"
    truncated = len(df) > max_rows
    display = df.head(max_rows)
    header = "".join(f"<th>{c}</th>" for c in display.columns)
    body = ""
    for _, row in display.iterrows():
        body += "<tr>" + "".join(f"<td>{v}</td>" for v in row) + "</tr>\n"
    note = f"<p><em>Showing first {max_rows} of {len(df)} flagged samples.</em></p>" if truncated else ""
    return note + f"""
    <table class="data-table">
      <thead><tr>{header}</tr></thead>
      <tbody>{body}</tbody>
    </table>"""

def concordance_table_html(df, flag_threshold, max_rows=200):
    """Render concordance pairs table with red highlights for duplicates."""
    if df is None or len(df) == 0:
        return ""
    truncated = len(df) > max_rows
    display = df.head(max_rows)
    header = "<tr><th>Sample A</th><th>Sample B</th><th>SNPs Called</th><th>SNPs Matched</th><th>Concordance %</th></tr>"
    body = ""
    for _, row in display.iterrows():
        pct = row["CONCORDANCE_PCT"]
        is_dup = pct >= flag_threshold
        row_style = ' style="background:#FFEBEE;"' if is_dup else ""
        flag_icon = " 🔴" if is_dup else ""
        body += (f"<tr{row_style}>"
                 f"<td>{row['SAMPLE_A']}</td>"
                 f"<td>{row['SAMPLE_B']}</td>"
                 f"<td>{int(row['N_CALLED']):,}</td>"
                 f"<td>{int(row['N_MATCH']):,}</td>"
                 f"<td><strong>{pct:.2f}%{flag_icon}</strong></td>"
                 f"</tr>\n")
    note = f"<p><em>Showing first {max_rows} of {len(df)} pairs.</em></p>" if truncated else ""
    return note + f"""
    <table class="data-table">
      <thead>{header}</thead>
      <tbody>{body}</tbody>
    </table>"""

def contamination_table_html(concordance_flag_pairs, samplesheet, flag_threshold):
    """
    Table of contaminated/duplicate samples (those in pairs >= flag_threshold).
    Shows sample ID, plate, well, and which sample(s) it's duplicated with.
    """
    if concordance_flag_pairs is None or len(concordance_flag_pairs) == 0:
        return "<p class='ok'>✓ No contaminated/duplicate samples detected (no pairs ≥ {:.0f}%).</p>".format(flag_threshold)

    # Build lookup: sample_id → plate, well
    lookup = {}
    if samplesheet is not None and has_plate_info:
        for _, srow in samplesheet.iterrows():
            sid = str(srow["sample_id"]).strip()
            lookup[sid] = {
                "plate": str(srow.get("plate", "N/A")),
                "well":  str(srow.get("well",  "N/A")),
            }

    # Collect all flagged samples and their duplicate partners
    dup_info = {}   # sample_id → {partners, max_conc, plate, well}
    for _, row in concordance_flag_pairs.iterrows():
        sa, sb, pct = str(row["SAMPLE_A"]), str(row["SAMPLE_B"]), row["CONCORDANCE_PCT"]
        for primary, partner in [(sa, sb), (sb, sa)]:
            if primary not in dup_info:
                dup_info[primary] = {"partners": [], "max_conc": 0.0,
                                     "plate": lookup.get(primary, {}).get("plate", "N/A"),
                                     "well":  lookup.get(primary, {}).get("well",  "N/A")}
            dup_info[primary]["partners"].append(f"{partner} ({pct:.2f}%)")
            dup_info[primary]["max_conc"] = max(dup_info[primary]["max_conc"], pct)

    rows_html = ""
    for sid, info in sorted(dup_info.items(), key=lambda x: -x[1]["max_conc"]):
        partners_str = "; ".join(info["partners"])
        rows_html += (f"<tr style='background:#FFEBEE;'>"
                      f"<td><strong>{sid}</strong></td>"
                      f"<td>{info['plate']}</td>"
                      f"<td>{info['well']}</td>"
                      f"<td>{info['max_conc']:.2f}%</td>"
                      f"<td>{partners_str}</td>"
                      f"</tr>\n")

    return f"""
    <table class="data-table">
      <thead><tr>
        <th>Sample ID</th><th>Plate</th><th>Well</th>
        <th>Max Concordance</th><th>Duplicate Partner(s)</th>
      </tr></thead>
      <tbody>{rows_html}</tbody>
    </table>"""

# ══════════════════════════════════════════════════════════════════════════════
# 6.  BUILD CONCORDANCE SECTION HTML
# ══════════════════════════════════════════════════════════════════════════════

def build_concordance_section():
    if concordance_df is None:
        return ""   # no concordance file provided, skip section entirely

    # --- Summary stats block ---
    total_pairs = len(concordance_df)
    section = f"""
<div class="card" id="concordance">
  <div class="card-header"><h2>🔗 Genotype Concordance – Sample Similarity</h2></div>
  <div class="card-body">
    <p>
      Total pairs evaluated: <strong>{total_pairs:,}</strong>
      &nbsp;|&nbsp; Pairs ≥ {args.concordance_warn}%: <strong>{n_warn_pairs}</strong>
      &nbsp;|&nbsp; Likely duplicates (≥ {args.concordance_flag}%): <strong style="color:#C62828">{n_flag_pairs}</strong>
    </p>
    <p style="margin-top:6px; font-size:13px; color:#546E7A;">
      <em>Unrelated samples from a shared-ancestry population typically share 50–65% genotypes.
      Pairs above {args.concordance_warn}% are shown below. Pairs ≥ {args.concordance_flag}%
      (highlighted in red 🔴) are likely duplicate samples or identical twins.</em>
    </p>
"""

    # --- Distribution plot ---
    if plot_concordance_b64:
        section += f'    <img class="plot-img" src="data:image/png;base64,{plot_concordance_b64}" alt="Concordance distribution" style="margin-top:14px;">\n'

    # --- Table or "no contamination" message ---
    if n_warn_pairs == 0:
        section += f"""
    <div style="margin-top:18px; padding:14px 18px; background:#E8F5E9; border-radius:6px; border:1px solid #A5D6A7;">
      <strong style="color:#2E7D32;">✓ No contamination detected.</strong>
      No sample pairs exceeded the {args.concordance_warn}% similarity threshold.
    </div>
"""
    else:
        section += f"""
    <h3 style="font-size:0.95rem; margin:18px 0 8px;">Pairs ≥ {args.concordance_warn}% Concordance</h3>
    {concordance_table_html(concordance_warn_pairs, args.concordance_flag)}
"""

    section += "  </div>\n</div>\n"
    return section


def build_plate_layout_section():
    if not plate_layout_b64s:
        return ""

    section = """
<div class="card" id="plate_layout">
  <div class="card-header"><h2>🧫 Plate Layouts – Concordance by Well</h2></div>
  <div class="card-body">
    <p style="font-size:13px; color:#546E7A; margin-bottom:14px;">
      Each well shows the sample ID and its highest pairwise concordance with any other sample.
      Wells highlighted in <span style="color:#B71C1C; font-weight:600;">red border</span>
      exceed the duplicate threshold (≥ {flag}%).<br>
      The <strong>coloured triangle</strong> in the top-right corner of each well indicates sex-check status:
      <span style="color:#4CAF50; font-weight:600;">▲ green</span> = sex OK,
      <span style="color:#F44336; font-weight:600;">▲ red</span> = sex DISCORDANT (collected ≠ inferred),
      <span style="color:#FF9800; font-weight:600;">▲ amber</span> = ambiguous F-statistic,
      <span style="color:#9E9E9E; font-weight:600;">▲ grey</span> = no sex-check data.
      Green wells (background) are below the {warn}% concordance reporting threshold.
    </p>
""".format(flag=args.concordance_flag, warn=args.concordance_warn)

    for plate_name, b64 in plate_layout_b64s.items():
        section += f"""
    <h3 style="font-size:0.95rem; margin:16px 0 6px; color:#283593;">{plate_name}</h3>
    <img class="plot-img" src="data:image/png;base64,{b64}" alt="Plate layout {plate_name}" style="max-width:100%;">
"""

    section += "  </div>\n</div>\n"
    return section


def build_plate_gc_section():
        if not plate_gc_b64s:
                return ""
        section = """
<div class="card" id="plate_gc">
    <div class="card-header"><h2>🧫 Plate QC – Call Rate vs 10% GenCall</h2></div>
    <div class="card-body">
        <p style="font-size:13px; color:#546E7A; margin-bottom:14px;">
            Per-plate scatter plots of Call Rate vs 10% GenCall score. Plates with
            poor hybridisation or failed processing will cluster in the lower-left
            corner (low call rate, low GenCall). Use these to spot plate-level failures.
        </p>
"""
        for plate_name, b64 in plate_gc_b64s.items():
                section += f"\n    <h3 style=\"font-size:0.95rem; margin:16px 0 6px; color:#283593;\">{plate_name}</h3>"
                section += f"\n    <img class=\"plot-img\" src=\"data:image/png;base64,{b64}\" alt=\"Plate GC {plate_name}\" style=\"max-width:100%;\">"

        section += "\n  </div>\n</div>\n"
        return section


def build_contamination_section():
    if concordance_flag_pairs is None:
        return ""

    section = f"""
<div class="card" id="contamination">
  <div class="card-header"><h2>⚠️ Contaminated / Duplicate Samples</h2></div>
  <div class="card-body">
    <p style="margin-bottom:12px; font-size:13px; color:#546E7A;">
      Samples listed here share ≥ {args.concordance_flag}% genotype concordance with at least one
      other sample in the cohort. This likely indicates sample duplication, contamination,
      or a labelling error. Both samples in each pair are listed.
    </p>
    {contamination_table_html(concordance_flag_pairs, samplesheet, args.concordance_flag)}
  </div>
</div>
"""
    return section

# ══════════════════════════════════════════════════════════════════════════════
# 7.  ASSEMBLE FINAL HTML
# ══════════════════════════════════════════════════════════════════════════════

now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

concordance_nav  = '<a href="#concordance">Concordance</a>' if concordance_df is not None else ""
plate_nav        = '<a href="#plate_layout">Plate Layouts</a>' if plate_layout_b64s else ""
contam_nav       = '<a href="#contamination">Contamination</a>' if (concordance_flag_pairs is not None and len(concordance_flag_pairs)) else ""
plate_gc_nav     = '<a href="#plate_gc">Plate QC (CallRate vs GenCall)</a>' if plate_gc_b64s else ""

# GTC QC explanatory note (shown when poor GC10 samples exist)
gtc_qc_note = ""
if n_poor_gc10:
    gtc_qc_note = (
        f"<p style=\"color:#B71C1C; font-weight:600; margin-top:8px;\">"
        f"GTC QC: <strong>{int(n_poor_gc10)}</strong> samples have low 10% GenCall (p10_gc). "
        "A low p10_gc means a meaningful fraction of that sample's SNP calls have poor cluster quality — "
        "this is often missed by PLINK call-rate QC alone, because PLINK counts any non-missing genotype as present "
        "even when the underlying cluster was poor. "
        "These samples are flagged as <strong>POOR_GC10</strong> in the flagged-samples table and should be inspected or excluded.</p>"
    )

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>H3aFlow – QC Report</title>
<style>
  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    font-size: 14px; background: #F5F7FA; color: #212121; line-height: 1.55;
  }}
  header {{
    background: linear-gradient(135deg, #1565C0 0%, #283593 100%);
    color: white; padding: 28px 40px 22px;
  }}
  header h1 {{ font-size: 1.7rem; font-weight: 700; letter-spacing: -0.3px; }}
  header p  {{ margin-top: 4px; opacity: 0.82; font-size: 0.9rem; }}
  main {{ max-width: 1200px; margin: 30px auto; padding: 0 24px 60px; }}
  .card {{
    background: white; border-radius: 8px;
    box-shadow: 0 1px 4px rgba(0,0,0,0.10); margin-bottom: 28px; overflow: hidden;
  }}
  .card-header {{ background: #E8EAF6; padding: 14px 22px; border-bottom: 1px solid #C5CAE9; }}
  .card-header h2 {{ font-size: 1.05rem; font-weight: 600; color: #283593; }}
  .card-body {{ padding: 20px 22px; }}
  .summary-table, .data-table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
  .summary-table th, .data-table th {{
    background: #283593; color: white; padding: 8px 12px; text-align: left; font-weight: 600;
  }}
  .summary-table td, .data-table td {{ padding: 7px 12px; border-bottom: 1px solid #EEEEEE; }}
  .summary-table tr:last-child td, .data-table tr:last-child td {{ border-bottom: none; }}
  .summary-table tr:nth-child(even) td,
  .data-table    tr:nth-child(even) td {{ background: #F8F9FE; }}
  .ok {{ color: #2E7D32; font-weight: 600; padding: 10px 0; }}
  .plot-img {{ max-width: 100%; border-radius: 6px; border: 1px solid #E0E0E0; margin-top: 10px; }}
  .params-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(220px, 1fr)); gap: 12px; margin-top: 6px; }}
  .param-box {{ background: #F3F4FB; border: 1px solid #C5CAE9; border-radius: 6px; padding: 10px 14px; }}
  .param-box .label {{ font-size: 11px; color: #5C6BC0; text-transform: uppercase; letter-spacing: 0.5px; }}
  .param-box .value {{ font-size: 1.1rem; font-weight: 700; color: #283593; margin-top: 2px; }}
  nav.toc {{ display: flex; flex-wrap: wrap; gap: 8px; margin: 0 0 24px; }}
  nav.toc a {{
    background: #E8EAF6; color: #283593; padding: 5px 14px; border-radius: 20px;
    text-decoration: none; font-size: 13px; font-weight: 500; transition: background 0.15s;
  }}
  nav.toc a:hover {{ background: #C5CAE9; }}
  footer {{ text-align: center; color: #9E9E9E; font-size: 12px; padding: 20px 0 40px; }}
</style>
</head>
<body>

<header>
  <h1>H3aFlow – QC Report</h1>
  <p>Generated: {now} &nbsp;|&nbsp; {total_samples} samples</p>
</header>

<main>

<nav class="toc">
  <a href="#overview">Overview</a>
  <a href="#sampleqc">Sample QC</a>
  <a href="#sexcheck">Sex Check</a>
  <a href="#ibd">IBD</a>
  <a href="#pca">PCA</a>
  {concordance_nav}
  {plate_nav}
  {plate_gc_nav}
  {contam_nav}
  <a href="#flagged">Flagged Samples</a>
</nav>

<div class="card" id="overview">
  <div class="card-header"><h2>📊 Cohort Overview</h2></div>
  <div class="card-body">
    <div class="params-grid">
      <div class="param-box"><div class="label">Total Samples</div><div class="value">{total_samples}</div></div>
      <div class="param-box"><div class="label">Samples Pass QC</div><div class="value">{int(n_pass_any)}</div></div>
      <div class="param-box"><div class="label">Uniquely Flagged</div><div class="value">{len(flagged_agg)}</div></div>
      <div class="param-box"><div class="label">Sex Discordant</div><div class="value">{n_sex_discord}</div></div>
      <div class="param-box"><div class="label">IBD Pairs Flagged</div><div class="value">{n_ibd_pairs}</div></div>
      <div class="param-box"><div class="label">PCA Outliers</div><div class="value">{n_pca_outliers}</div></div>
      <div class="param-box"><div class="label">Concordance ≥{args.concordance_warn}%</div><div class="value">{n_warn_pairs}</div></div>
      <div class="param-box"><div class="label" style="color:#B71C1C;">Likely Duplicates (≥{args.concordance_flag}%)</div><div class="value" style="color:#B71C1C;">{n_flag_pairs}</div></div>
    </div>
    <br>
    <h3 style="font-size:0.95rem; margin-bottom:8px;">QC Parameters Used</h3>
    <div class="params-grid">
      <div class="param-box"><div class="label">Missingness cutoff (--mind)</div><div class="value">{args.mind}</div></div>
      <div class="param-box"><div class="label">Het SD cutoff</div><div class="value">±{args.het_sd} SD</div></div>
      <div class="param-box"><div class="label">IBD PI_HAT threshold</div><div class="value">{args.pi_hat}</div></div>
      <div class="param-box"><div class="label">Concordance report threshold</div><div class="value">{args.concordance_warn}%</div></div>
      <div class="param-box"><div class="label">Concordance duplicate flag</div><div class="value">{args.concordance_flag}%</div></div>
    </div>
    <br>
    {summary_table_html(summary_df)}
  </div>
</div>

<div class="card" id="sampleqc">
  <div class="card-header"><h2>🔬 Sample QC – Call Rate &amp; Heterozygosity</h2></div>
  <div class="card-body">
    <p>Samples failing call-rate: <strong>{int(n_fail_mind)}</strong>
       &nbsp;|&nbsp; Failing heterozygosity: <strong>{int(n_fail_het)}</strong>
       &nbsp;|&nbsp; Failing either: <strong>{int(n_fail_any)}</strong></p>
    <img class="plot-img" src="data:image/png;base64,{plot_sampleqc_b64}" alt="Sample QC plots">
        {gtc_qc_note}
  </div>
</div>

<div class="card" id="sexcheck">
  <div class="card-header"><h2>⚥ Sex Check – {"GTC computed gender" if use_gtc_sex else "chrX F-Statistic"} &amp; XY Intensities</h2></div>
  <div class="card-body">
    <p>PLINK STATUS=PROBLEM: <strong>{n_sex_problem}</strong>
       &nbsp;|&nbsp; Discordant (inferred ≠ collected): <strong>{n_sex_discord}</strong></p>
    <img class="plot-img" src="data:image/png;base64,{plot_sexcheck_b64}" alt="Sex check plots">
  </div>
</div>

<div class="card" id="ibd">
  <div class="card-header"><h2>🧬 IBD – Identity by Descent / Relatedness</h2></div>
  <div class="card-body">
    <p>Pairs with PI_HAT ≥ {args.pi_hat}: <strong>{n_ibd_pairs}</strong>
       &nbsp;|&nbsp; Duplicates / MZ twins: <strong>{n_ibd_dup}</strong></p>
    <img class="plot-img" src="data:image/png;base64,{plot_ibd_b64}" alt="IBD plots">
  </div>
</div>

<div class="card" id="pca">
  <div class="card-header"><h2>📐 PCA – Population Structure</h2></div>
  <div class="card-body">
    <p>PCA outliers (&gt;6 SD on PC1/PC2): <strong>{n_pca_outliers}</strong></p>
    <img class="plot-img" src="data:image/png;base64,{plot_pca_b64}" alt="PCA plots">
  </div>
</div>

{build_concordance_section()}
{build_plate_layout_section()}
{build_plate_gc_section()}
{build_contamination_section()}

<div class="card" id="flagged">
  <div class="card-header"><h2>🚩 Flagged Samples (union of all QC stages)</h2></div>
  <div class="card-body">
    {flagged_table_html(flagged_agg)}
  </div>
</div>

</main>

<footer>H3aFlow QC Report &nbsp;·&nbsp; {now}</footer>

</body>
</html>"""

with open(args.out_html, "w") as fh:
    fh.write(html)

print(f"[generate_report] HTML report → {args.out_html}", flush=True)
print(f"[generate_report] Done. {len(flagged_agg)} samples flagged across all QC stages.")