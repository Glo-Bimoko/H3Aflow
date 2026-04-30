"""
plot_pca.py
===========
Reads PLINK2 --pca output files (.eigenvec and .eigenval) and produces a
single PNG with:
  Panel 1: PC1 vs PC2 scatter (coloured by density; outliers annotated)
  Panel 2: PC1 vs PC3 scatter
  Panel 3: Scree plot – % variance explained per PC

Usage (called from pca.nf):
  python plot_pca.py \
      --eigenvec pca.eigenvec \
      --eigenval pca.eigenval \
      --out      pca_plot.png

Optional:
  --sex_info  path/to/sex_info.tsv   (adds sex-colouring if provided)
  --n_pcs     int                    (number of PCs to show in scree; default: all)
  --pc_x      int                    (PC for x-axis of scatter; default: 1)
  --pc_y      int                    (PC for y-axis of scatter; default: 2)
  --pc_z      int                    (PC for second scatter y-axis; default: 3)
  --sd_cutoff float                  (flag outliers beyond N SD on PC1+PC2; default: 6)

PLINK2 .eigenvec format:
  #FID IID PC1 PC2 … PCn
  (header line starts with '#')
"""

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize

# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Plot PLINK2 PCA results")
parser.add_argument("--eigenvec",  required=True)
parser.add_argument("--eigenval",  required=True)
parser.add_argument("--out",       required=True)
parser.add_argument("--sex_info",  default=None,
                    help="Optional TSV: sampleid, sex (0=Female, 1=Male)")
parser.add_argument("--n_pcs",     type=int, default=None,
                    help="Number of PCs to display in scree plot")
parser.add_argument("--pc_x",      type=int, default=1)
parser.add_argument("--pc_y",      type=int, default=2)
parser.add_argument("--pc_z",      type=int, default=3)
parser.add_argument("--sd_cutoff", type=float, default=6.0,
                    help="Flag samples beyond N SD on PC1 and PC2 (default: 6)")
args = parser.parse_args()

# ── Load eigenvec ──────────────────────────────────────────────────────────────
print(f"[plot_pca] Reading {args.eigenvec} …", flush=True)
evec = pd.read_csv(args.eigenvec, sep=r"\s+", comment=None)

# PLINK2 adds '#FID' as first column header
evec.columns = evec.columns.str.lstrip("#")
if "FID" not in evec.columns and "IID" not in evec.columns:
    sys.exit("[plot_pca] ERROR: .eigenvec has unexpected column structure.")

# Identify PC columns
pc_cols = [c for c in evec.columns if c.upper().startswith("PC")]
if not pc_cols:
    sys.exit("[plot_pca] ERROR: No PC columns found in .eigenvec.")

n_pcs_available = len(pc_cols)
print(f"[plot_pca] {len(evec)} samples, {n_pcs_available} PCs", flush=True)

def pc_col(n):
    """Return column name for PC n (1-indexed)."""
    name = f"PC{n}"
    if name not in evec.columns:
        sys.exit(f"[plot_pca] ERROR: PC{n} not found in .eigenvec.")
    return name

# ── Load eigenval ──────────────────────────────────────────────────────────────
print(f"[plot_pca] Reading {args.eigenval} …", flush=True)
evals = pd.read_csv(args.eigenval, header=None, names=["eigenvalue"])
evals["pct_var"] = 100.0 * evals["eigenvalue"] / evals["eigenvalue"].sum()
evals["cumvar"]  = evals["pct_var"].cumsum()
n_scree = args.n_pcs if args.n_pcs else len(evals)
evals_plot = evals.head(n_scree).copy()
evals_plot["PC"] = [f"PC{i+1}" for i in range(len(evals_plot))]

# ── Optional sex colouring ─────────────────────────────────────────────────────
SEX_PALETTE = {"Male": "#2196F3", "Female": "#E91E63", "Unknown": "#9E9E9E"}
use_sex = False
if args.sex_info:
    try:
        sex_info = pd.read_csv(args.sex_info, sep="\t", dtype=str)
        sex_info.columns = sex_info.columns.str.strip().str.lower()
        sex_info = sex_info.rename(columns={"sampleid": "IID"})
        sex_info["IID"] = sex_info["IID"].str.strip()
        sex_info["COLLECTED_SEX"] = sex_info["sex"].map(
            {"0": "Female", "1": "Male"}
        ).fillna("Unknown")
        evec = evec.merge(sex_info[["IID", "COLLECTED_SEX"]], on="IID", how="left")
        evec["COLLECTED_SEX"] = evec["COLLECTED_SEX"].fillna("Unknown")
        use_sex = True
        print("[plot_pca] Sex information loaded for colouring.", flush=True)
    except Exception as e:
        print(f"[plot_pca] WARNING: Could not load sex_info ({e}); using density colouring.", flush=True)

# ── Identify PC outliers ───────────────────────────────────────────────────────
col_x = pc_col(args.pc_x)
col_y = pc_col(args.pc_y)

for col in [col_x, col_y]:
    mean_ = evec[col].mean()
    sd_   = evec[col].std()
    evec[f"{col}_outlier"] = (
        (evec[col] > mean_ + args.sd_cutoff * sd_) |
        (evec[col] < mean_ - args.sd_cutoff * sd_)
    )

evec["PCA_OUTLIER"] = evec[f"{col_x}_outlier"] | evec[f"{col_y}_outlier"]
n_outliers = evec["PCA_OUTLIER"].sum()
print(f"[plot_pca] PCA outliers (>{args.sd_cutoff} SD on {col_x} or {col_y}): {n_outliers}", flush=True)

# ── Helper: scatter panel ──────────────────────────────────────────────────────
def scatter_panel(ax, col_a, col_b, title, eigenval_df):
    pct_a = eigenval_df.loc[eigenval_df["PC"] == col_a, "pct_var"].values
    pct_b = eigenval_df.loc[eigenval_df["PC"] == col_b, "pct_var"].values
    xlabel = f"{col_a}  ({pct_a[0]:.2f}% var)" if len(pct_a) else col_a
    ylabel = f"{col_b}  ({pct_b[0]:.2f}% var)" if len(pct_b) else col_b

    main = evec[~evec["PCA_OUTLIER"]]
    outs = evec[evec["PCA_OUTLIER"]]

    if use_sex:
        for sex, grp in main.groupby("COLLECTED_SEX"):
            ax.scatter(grp[col_a], grp[col_b],
                       c=SEX_PALETTE.get(sex, "#9E9E9E"),
                       alpha=0.55, s=12, edgecolors="none", label=sex)
        if len(outs):
            ax.scatter(outs[col_a], outs[col_b],
                       c="black", marker="x", s=40, linewidths=0.8,
                       label=f"Outlier (>{args.sd_cutoff}σ)")
        ax.legend(title="Collected Sex", fontsize=7, framealpha=0.7)
    else:
        # Density colouring via histogram2d
        x = main[col_a].values
        y = main[col_b].values
        if len(x) > 1:
            h, xedges, yedges = np.histogram2d(x, y, bins=100)
            xidx = np.clip(np.searchsorted(xedges, x) - 1, 0, h.shape[0]-1)
            yidx = np.clip(np.searchsorted(yedges, y) - 1, 0, h.shape[1]-1)
            density = h[xidx, yidx]
            order = density.argsort()
            ax.scatter(x[order], y[order], c=density[order],
                       cmap="viridis", alpha=0.6, s=12, edgecolors="none")
        else:
            ax.scatter(x, y, c="#1976D2", alpha=0.6, s=12)
        if len(outs):
            ax.scatter(outs[col_a], outs[col_b],
                       c="red", marker="x", s=40, linewidths=0.8,
                       label=f"Outlier (>{args.sd_cutoff}σ)")
            ax.legend(fontsize=7, framealpha=0.7)

    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.grid(True, linewidth=0.4, alpha=0.5)

# ── Build figure ───────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle("Principal Component Analysis (PCA) – QC summary", fontsize=14, fontweight="bold")

# Panel 1: PC_x vs PC_y
scatter_panel(axes[0], col_x, col_y,
              f"{col_x} vs {col_y}", evals_plot)

# Panel 2: PC_x vs PC_z (only if PC_z exists)
col_z = f"PC{args.pc_z}"
if col_z in evec.columns:
    scatter_panel(axes[1], col_x, col_z,
                  f"{col_x} vs {col_z}", evals_plot)
else:
    axes[1].set_visible(False)

# Panel 3: Scree plot
ax = axes[2]
bar_colors = plt.cm.Blues(np.linspace(0.4, 0.85, len(evals_plot)))
bars = ax.bar(evals_plot["PC"], evals_plot["pct_var"],
              color=bar_colors, edgecolor="white", linewidth=0.5)
ax2 = ax.twinx()
ax2.plot(evals_plot["PC"], evals_plot["cumvar"],
         color="#E65100", marker="o", markersize=4, linewidth=1.5,
         label="Cumulative %")
ax2.set_ylabel("Cumulative % variance", fontsize=10, color="#E65100")
ax2.tick_params(axis="y", labelcolor="#E65100")
ax2.set_ylim(0, 105)
ax.set_xlabel("Principal Component", fontsize=10)
ax.set_ylabel("% Variance Explained", fontsize=10)
ax.set_title("Scree Plot", fontsize=11)
ax.tick_params(axis="x", rotation=45)
ax.grid(axis="y", linewidth=0.4, alpha=0.5)

plt.tight_layout()
plt.savefig(args.plot, dpi=150, bbox_inches="tight")
plt.close()
print(f"[plot_pca] Plot saved → {args.plot}", flush=True)

# ── Print variance table ───────────────────────────────────────────────────────
print("\n[plot_pca] Variance explained:")
for _, row in evals_plot.iterrows():
    print(f"  {row['PC']:6s}: {row['pct_var']:6.2f}%  (cumulative: {row['cumvar']:6.2f}%)")

if n_outliers:
    print(f"\n[plot_pca] PCA outlier IIDs (>{args.sd_cutoff} SD):")
    for iid in evec.loc[evec["PCA_OUTLIER"], "IID"].tolist():
        print(f"  {iid}")
