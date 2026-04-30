"""
flag_ibd_duplicates.py
======================
Reads the PLINK --genome output (ibd.genome), flags pairs exceeding the
PI_HAT threshold, and writes:
  --out   ibd_duplicates.tsv   (one row per flagged pair + relationship label)
  --plot  ibd_plot.png         (Z0 vs Z1 scatter; PI_HAT histogram)

Relationship labels (approximate, based on IBD sharing):
  PI_HAT ≥ 0.9  → Duplicate / MZ twin
  PI_HAT ≥ 0.45 → 1st degree (parent-offspring or full sibling)
  PI_HAT ≥ 0.20 → 2nd degree (half-sibling, grandparent, aunt/uncle)
  below cutoff   → Not flagged (not written to output)

The --min flag passed to PLINK already restricts the .genome file to pairs
above --pi_hat, so we re-apply the threshold here for safety and labelling.

PLINK .genome columns of interest:
  FID1 IID1 FID2 IID2 RT EZ Z0 Z1 Z2 PI_HAT PHE DST PPC RATIO
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Flag IBD duplicates from PLINK --genome")
parser.add_argument("--genome",  required=True,  help="PLINK .genome file")
parser.add_argument("--pi_hat",  type=float, default=0.1875,
                    help="PI_HAT threshold for flagging (default: 0.1875)")
parser.add_argument("--out",     required=True,  help="Output TSV path")
parser.add_argument("--plot",    required=True,  help="Output PNG path")
args = parser.parse_args()

# ── Load .genome ───────────────────────────────────────────────────────────────
print(f"[flag_ibd_duplicates] Reading {args.genome} …", flush=True)
genome = pd.read_csv(args.genome, sep=r"\s+")
genome.columns = genome.columns.str.strip()

required_cols = {"IID1", "IID2", "PI_HAT"}
if not required_cols.issubset(genome.columns):
    import sys
    sys.exit(
        f"[flag_ibd_duplicates] ERROR: Missing columns in .genome file.\n"
        f"  Required: {required_cols}\n  Found: {list(genome.columns)}"
    )

print(f"[flag_ibd_duplicates] Total pairs in .genome: {len(genome)}", flush=True)

# ── Filter to pairs above threshold ───────────────────────────────────────────
flagged = genome[genome["PI_HAT"] >= args.pi_hat].copy()
print(
    f"[flag_ibd_duplicates] Pairs with PI_HAT ≥ {args.pi_hat}: {len(flagged)}",
    flush=True,
)

# ── Label relationship ─────────────────────────────────────────────────────────
def relationship_label(pi):
    if pi >= 0.90:
        return "Duplicate/MZ_twin"
    elif pi >= 0.45:
        return "1st_degree"
    elif pi >= 0.20:
        return "2nd_degree"
    else:
        return "Cryptic_relatedness"

flagged["RELATIONSHIP"] = flagged["PI_HAT"].apply(relationship_label)

# ── Determine which sample to recommend removing ───────────────────────────────
# Heuristic: within each pair, flag IID2 for removal (keeps IID1 by default).
# Users can override; this is just a starting recommendation.
flagged["RECOMMEND_REMOVE"] = flagged["IID2"]

# Select and order output columns
out_cols = [c for c in
            ["FID1","IID1","FID2","IID2","Z0","Z1","Z2","PI_HAT","RELATIONSHIP","RECOMMEND_REMOVE"]
            if c in flagged.columns]
flagged[out_cols].to_csv(args.out, sep="\t", index=False)
print(f"[flag_ibd_duplicates] Written → {args.out}", flush=True)

# ── Relationship breakdown ─────────────────────────────────────────────────────
for rel, grp in flagged.groupby("RELATIONSHIP"):
    print(f"  {rel:25s}: {len(grp)} pairs")

# ── Plot ───────────────────────────────────────────────────────────────────────
REL_PALETTE = {
    "Duplicate/MZ_twin":    "#D32F2F",
    "1st_degree":           "#FF9800",
    "2nd_degree":           "#1976D2",
    "Cryptic_relatedness":  "#616161",
}

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle("IBD / Relatedness QC", fontsize=14, fontweight="bold")

# --- Panel 1: Z0 vs Z1 (IBD state proportions) --------------------------------
ax = axes[0]
# Plot all pairs in genome (below threshold) as grey background
if len(genome) > len(flagged):
    background = genome[genome["PI_HAT"] < args.pi_hat]
    if "Z0" in background.columns and "Z1" in background.columns:
        ax.scatter(
            background["Z0"], background["Z1"],
            c="#BDBDBD", alpha=0.3, s=8, edgecolors="none", label="Below threshold"
        )

# Overlay flagged pairs coloured by relationship
for rel, grp in flagged.groupby("RELATIONSHIP"):
    if "Z0" in grp.columns and "Z1" in grp.columns:
        ax.scatter(
            grp["Z0"], grp["Z1"],
            c=REL_PALETTE.get(rel, "#9E9E9E"),
            alpha=0.7, s=30, edgecolors="none", label=rel
        )

ax.set_xlabel("Z0  (P(IBD=0))", fontsize=11)
ax.set_ylabel("Z1  (P(IBD=1))", fontsize=11)
ax.set_title(f"IBD State Proportions\n(PI_HAT threshold = {args.pi_hat})", fontsize=11)
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-0.05, 1.05)
ax.legend(fontsize=8, framealpha=0.7)
ax.grid(True, linewidth=0.4, alpha=0.5)

# --- Panel 2: PI_HAT histogram (flagged pairs only) ---------------------------
ax = axes[1]
bins = np.linspace(args.pi_hat, 1.0, 30)
for rel, grp in flagged.groupby("RELATIONSHIP"):
    ax.hist(
        grp["PI_HAT"],
        bins=bins,
        color=REL_PALETTE.get(rel, "#9E9E9E"),
        alpha=0.75,
        label=rel,
        edgecolor="white",
        linewidth=0.4,
    )

ax.axvline(args.pi_hat, color="black", linestyle="--", linewidth=1.2,
           label=f"Threshold ({args.pi_hat})")
ax.set_xlabel("PI_HAT", fontsize=11)
ax.set_ylabel("Number of pairs", fontsize=11)
ax.set_title("PI_HAT Distribution\n(flagged pairs)", fontsize=11)
ax.legend(fontsize=8, framealpha=0.7)
ax.grid(True, linewidth=0.4, alpha=0.5)

plt.tight_layout()
plt.savefig(args.plot, dpi=150, bbox_inches="tight")
plt.close()
print(f"[flag_ibd_duplicates] Plot saved → {args.plot}", flush=True)
