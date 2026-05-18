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
  --sexcheck    sex_check.sexcheck          (PLINK --check-sex; annotated by annotate_sex_check.py)
  --xy_tsv      xy_intensity.tsv           (extract_xy_intensity.py)
  --qc_stats    sample_qc_stats.tsv        (compute_sample_qc.py)
  --genome      ibd.genome                 (PLINK --genome)
  --eigenvec    pca.eigenvec               (PLINK2 --pca)
  --sex_info    sex_info.tsv               (raw collected sex)
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

# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Generate HTML QC report")
parser.add_argument("--sexcheck",    required=True)
parser.add_argument("--xy_tsv",     required=True)
parser.add_argument("--qc_stats",   required=True)
parser.add_argument("--genome",     required=True)
parser.add_argument("--eigenvec",   required=True)
parser.add_argument("--sex_info",   required=True)
parser.add_argument("--out_html",   required=True)
parser.add_argument("--out_flagged",required=True)
parser.add_argument("--out_summary",required=True)
# Optional thresholds (match nextflow.config defaults)
parser.add_argument("--pi_hat",     type=float, default=0.1875)
parser.add_argument("--mind",       type=float, default=0.05)
parser.add_argument("--het_sd",     type=float, default=3.0)
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

# annotate_sex_check.py already adds COLLECTED_SEX, INFERRED_SEX, DISCORDANT.
# Only fall back to merging from sex_info if the column is genuinely absent
# (e.g. raw PLINK output was passed instead of the annotated file).
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

# PCA outliers (>6 SD on PC1/PC2)
for col in ["PC1","PC2"]:
    if col in evec.columns:
        m, s = evec[col].mean(), evec[col].std()
        evec[f"{col}_out"] = (evec[col].abs() > m + 6*s) | (evec[col] < m - 6*s)
evec["PCA_OUTLIER"] = (
    evec.get("PC1_out", pd.Series(False, index=evec.index)) |
    evec.get("PC2_out", pd.Series(False, index=evec.index))
)
n_pca_outliers = int(evec["PCA_OUTLIER"].sum())

# ══════════════════════════════════════════════════════════════════════════════
# 2.  BUILD FLAGGED SAMPLES TABLE
# ══════════════════════════════════════════════════════════════════════════════

print("[generate_report] Building flagged samples list …", flush=True)

flag_records = []

# Sample QC failures
for _, row in qc[qc["FAIL_ANY"]].iterrows():
    flag_records.append({"IID": str(row["IID"]), "FLAG": row["FAIL_REASON"], "SOURCE": "SAMPLE_QC"})

# Sex discordance
if "DISCORDANT" in sexcheck.columns:
    for _, row in sexcheck[sexcheck["DISCORDANT"]].iterrows():
        flag_records.append({"IID": str(row["IID"]), "FLAG": "SEX_DISCORDANT", "SOURCE": "SEX_CHECK"})

# IBD duplicates — flag IID2 in each pair
if len(flagged_ibd) and "IID2" in flagged_ibd.columns:
    for _, row in flagged_ibd.iterrows():
        flag_records.append({
            "IID": str(row["IID2"]),
            "FLAG": f"IBD_{row['RELATIONSHIP']}",
            "SOURCE": "IBD",
        })

# PCA outliers
for iid in evec.loc[evec["PCA_OUTLIER"], "IID"].tolist():
    flag_records.append({"IID": str(iid), "FLAG": "PCA_OUTLIER", "SOURCE": "PCA"})

flagged_df = pd.DataFrame(flag_records) if flag_records else pd.DataFrame(columns=["IID","FLAG","SOURCE"])

# Collapse multiple flags per sample
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
    ("Total samples",              total_samples),
    ("Samples passing call-rate QC", int(n_pass_mind)),
    ("Samples failing call-rate QC", int(n_fail_mind)),
    ("Samples failing heterozygosity QC", int(n_fail_het)),
    ("Samples failing ANY QC",     int(n_fail_any)),
    ("Sex-discordant samples",     n_sex_discord),
    ("PLINK sex STATUS=PROBLEM",   n_sex_problem),
    ("IBD pairs flagged",          n_ibd_pairs),
    ("  of which duplicates/MZ",   n_ibd_dup),
    ("PCA outliers (>6 SD)",       n_pca_outliers),
    ("Total uniquely flagged samples", len(flagged_agg)),
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

# ── Plot B: Sex check ─────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
fig.suptitle("Sex Check – chrX F-Statistic", fontsize=12, fontweight="bold")

ax = axes[0]
# COLLECTED_SEX is already on sexcheck (from annotate_sex_check.py).
# Use it directly — no merge needed.
for sex, grp in sexcheck.groupby("COLLECTED_SEX"):
    ax.hist(grp["F"].dropna(), bins=50, alpha=0.65,
            color=SEX_PALETTE.get(sex, "#9E9E9E"), label=sex,
            edgecolor="white", linewidth=0.3)
ax.axvline(0.2, color="grey", linestyle=":", linewidth=1)
ax.axvline(0.8, color="grey", linestyle=":", linewidth=1)
ax.set_xlabel("F-statistic (chrX inbreeding coefficient)", fontsize=10)
ax.set_ylabel("Samples", fontsize=10)
ax.set_title("F-statistic by collected sex", fontsize=10)
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

# ══════════════════════════════════════════════════════════════════════════════
# 5.  BUILD HTML
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

now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

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
    </div>
    <br>
    <h3 style="font-size:0.95rem; margin-bottom:8px;">QC Parameters Used</h3>
    <div class="params-grid">
      <div class="param-box"><div class="label">Missingness cutoff (--mind)</div><div class="value">{args.mind}</div></div>
      <div class="param-box"><div class="label">Het SD cutoff</div><div class="value">±{args.het_sd} SD</div></div>
      <div class="param-box"><div class="label">IBD PI_HAT threshold</div><div class="value">{args.pi_hat}</div></div>
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
  </div>
</div>

<div class="card" id="sexcheck">
  <div class="card-header"><h2>⚥ Sex Check – chrX F-Statistic &amp; XY Intensities</h2></div>
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