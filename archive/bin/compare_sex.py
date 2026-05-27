#!/usr/bin/env python3
"""
compare_sex.py
Compares PLINK --check-sex F-statistics against collected sex.

Changes from previous version (H3AGWAS alignment):
  - F thresholds configurable via --f_lo_male / --f_hi_female
  - Plot 1 changed from histogram to scatter: sample index vs F, coloured by
    collected sex. This is far more informative — you can immediately see
    whether the male/female clusters are well-separated and spot outliers.
    H3AGWAS visualises F distributions per sex group for the same reason.
  - Plot 2 (concordance bar chart) retained and improved.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

parser = argparse.ArgumentParser()
parser.add_argument("--plink_sex",     required=True)
parser.add_argument("--collected_sex", required=True)
parser.add_argument("--out",           required=True)
parser.add_argument("--plot",          required=True)
parser.add_argument("--f_lo_male",     type=float, default=0.8,
                    help="F lower bound for Male call (default: 0.8)")
parser.add_argument("--f_hi_female",   type=float, default=0.2,
                    help="F upper bound for Female call (default: 0.2)")
args = parser.parse_args()

f_lo_male   = args.f_lo_male
f_hi_female = args.f_hi_female

# ── Load PLINK sex check ──────────────────────────────────────────────────────
plink_sex = pd.read_csv(args.plink_sex, sep=r'\s+')
plink_sex["IID"] = plink_sex["IID"].astype(str).str.strip()

print(f"Loaded PLINK sex check with {len(plink_sex)} samples")
print(f"Columns: {plink_sex.columns.tolist()}")
print(f"STATUS value counts: {plink_sex['STATUS'].value_counts().to_dict() if 'STATUS' in plink_sex.columns else 'No STATUS column'}")
print(f"F thresholds — Male > {f_lo_male}, Female < {f_hi_female}")

if 'F' not in plink_sex.columns:
    print("Warning: No F-statistic column found in PLINK output")
    plink_sex['F'] = np.nan

# ── Load collected sex ────────────────────────────────────────────────────────
try:
    collected = pd.read_csv(args.collected_sex, sep="\t", dtype=str)
    collected.columns = [c.lower() for c in collected.columns]
    collected = collected.rename(columns={'sampleid': 'IID', 'sex': 'COLLECTED_SEX_NUM'})
    collected["IID"] = collected["IID"].astype(str).str.strip()
    collected['COLLECTED_SEX'] = collected['COLLECTED_SEX_NUM'].map(
        {'0': 'Female', '1': 'Male'}
    ).fillna('Unknown')
    print(f"Loaded collected sex for {len(collected)} samples")
    has_collected = True
except Exception as e:
    print(f"Warning: Could not read collected sex file: {e}")
    has_collected = False
    collected = pd.DataFrame({'IID': plink_sex['IID'], 'COLLECTED_SEX': 'Unknown'})

# ── Merge ─────────────────────────────────────────────────────────────────────
merged = plink_sex.merge(collected[['IID', 'COLLECTED_SEX']], on='IID', how='left')
merged['COLLECTED_SEX'] = merged['COLLECTED_SEX'].fillna('Unknown')

# ── Infer sex from F using configurable thresholds ────────────────────────────
if merged['F'].isna().all():
    print("All F-statistics are NaN — data may be monomorphic.")
    merged['INFERRED_SEX'] = 'Unknown'
else:
    def infer_sex(f_val):
        if pd.isna(f_val):     return 'Unknown'
        if f_val > f_lo_male:  return 'Male'
        if f_val < f_hi_female: return 'Female'
        return 'Ambiguous'
    merged['INFERRED_SEX'] = merged['F'].apply(infer_sex)

# ── Concordance ───────────────────────────────────────────────────────────────
def concordant(row):
    if row['INFERRED_SEX'] in ('Unknown', 'Ambiguous'):
        return 'Ambiguous/Unknown'
    if row['COLLECTED_SEX'] == 'Unknown':
        return 'No collected sex'
    return 'Concordant' if row['INFERRED_SEX'] == row['COLLECTED_SEX'] else 'Discordant'

merged['CONCORDANCE'] = merged.apply(concordant, axis=1)
merged['DISCORDANCE_TYPE'] = np.where(
    merged['CONCORDANCE'] == 'Discordant',
    merged['INFERRED_SEX'] + ' vs ' + merged['COLLECTED_SEX'],
    merged['CONCORDANCE']
)

# ── Save outputs ──────────────────────────────────────────────────────────────
discordant = merged[merged['CONCORDANCE'] == 'Discordant']
discordant.to_csv(args.out, sep="\t", index=False)
merged.to_csv(args.out.replace('.tsv', '_all.tsv').replace('.txt', '_all.tsv'),
              sep="\t", index=False)

# ── Plot ──────────────────────────────────────────────────────────────────────
# Colour scheme: collected sex determines point colour.
# Shape/edge distinguishes concordant vs discordant inference.
SEX_COLORS = {
    'Male':    '#2196F3',   # blue
    'Female':  '#E91E63',   # pink
    'Unknown': '#9E9E9E',   # grey
}

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# ── Panel 1: F-statistic scatter coloured by COLLECTED sex ───────────────────
# (H3AGWAS-style: show per-sex F distribution so cluster separation is visible)
ax1 = axes[0]

if merged['F'].notna().any():
    # Sort by collected sex for a clean visual grouping on x-axis
    order = {'Male': 0, 'Female': 1, 'Unknown': 2}
    plot_df = merged.copy()
    plot_df['_sex_order'] = plot_df['COLLECTED_SEX'].map(order).fillna(3)
    plot_df = plot_df.sort_values(['_sex_order', 'F']).reset_index(drop=True)

    for sex, color in SEX_COLORS.items():
        sub = plot_df[plot_df['COLLECTED_SEX'] == sex]
        if sub.empty:
            continue
        # Concordant → filled circle; discordant → red-edged diamond
        concord_mask = sub['CONCORDANCE'] == 'Concordant'

        ax1.scatter(sub[concord_mask].index, sub[concord_mask]['F'],
                    c=color, s=22, alpha=0.7, marker='o',
                    linewidths=0, label=f'Collected {sex}')
        discord_sub = sub[~concord_mask & sub['F'].notna()]
        if not discord_sub.empty:
            ax1.scatter(discord_sub.index, discord_sub['F'],
                        c=color, s=50, alpha=0.9, marker='D',
                        edgecolors='red', linewidths=1.2,
                        label=f'Discordant {sex}')

    # Threshold lines
    ax1.axhline(y=f_lo_male,   color='#1565C0', linestyle='--', linewidth=1.5,
                label=f'Male threshold ({f_lo_male})')
    ax1.axhline(y=f_hi_female, color='#AD1457', linestyle='--', linewidth=1.5,
                label=f'Female threshold ({f_hi_female})')
    ax1.axhspan(f_hi_female, f_lo_male, alpha=0.06, color='orange',
                label='Ambiguous zone')

    ax1.set_xlabel('Sample (sorted by collected sex)', fontsize=11)
    ax1.set_ylabel('F-statistic (chrX inbreeding coefficient)', fontsize=11)
    ax1.set_title('ChrX F-statistic by Collected Sex', fontsize=12, fontweight='bold')
    ax1.set_ylim(-0.1, 1.1)
    ax1.legend(fontsize=8, loc='upper right', ncol=1)
    ax1.grid(True, alpha=0.25)

    # Annotate sample counts per sex
    for sex, color in SEX_COLORS.items():
        n = (merged['COLLECTED_SEX'] == sex).sum()
        if n > 0:
            ax1.annotate(f'n={n}', xy=(0, 0), xycoords='axes fraction',
                         color=color, fontsize=8)
else:
    ax1.text(0.5, 0.5,
             'No F-statistics available\n(dataset may be monomorphic)',
             ha='center', va='center', transform=ax1.transAxes, fontsize=12)
    ax1.set_title('F-statistic — No Data', fontsize=12, fontweight='bold')

# ── Panel 2: Concordance summary bar chart ────────────────────────────────────
ax2 = axes[1]
if len(merged) > 0:
    cat_colors = {
        'Concordant':         '#43A047',
        'Discordant':         '#E53935',
        'Ambiguous/Unknown':  '#FB8C00',
        'No collected sex':   '#9E9E9E',
        'Male vs Female':     '#E53935',
        'Female vs Male':     '#E53935',
    }
    concordance_counts = merged['DISCORDANCE_TYPE'].value_counts()
    bar_colors = [cat_colors.get(x, '#757575') for x in concordance_counts.index]

    bars = ax2.bar(range(len(concordance_counts)), concordance_counts.values,
                   color=bar_colors, alpha=0.8, edgecolor='black', linewidth=0.7)
    ax2.set_xticks(range(len(concordance_counts)))
    ax2.set_xticklabels(concordance_counts.index, rotation=40, ha='right', fontsize=9)
    ax2.set_ylabel('Count', fontsize=11)
    ax2.set_title('Sex Concordance Summary', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')

    for bar, value in zip(bars, concordance_counts.values):
        ax2.text(bar.get_x() + bar.get_width() / 2., bar.get_height() + 0.3,
                 str(value), ha='center', va='bottom', fontsize=10, fontweight='bold')

    # Summary text box
    n_concord = (merged['CONCORDANCE'] == 'Concordant').sum()
    n_discord = (merged['CONCORDANCE'] == 'Discordant').sum()
    n_total   = len(merged)
    pct_c = 100 * n_concord / n_total if n_total > 0 else 0
    pct_d = 100 * n_discord / n_total if n_total > 0 else 0
    summary_txt = (f"Total: {n_total}\n"
                   f"Concordant: {n_concord} ({pct_c:.1f}%)\n"
                   f"Discordant: {n_discord} ({pct_d:.1f}%)")
    ax2.text(0.97, 0.97, summary_txt,
             transform=ax2.transAxes, fontsize=9,
             verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow',
                       edgecolor='grey', alpha=0.8))
else:
    ax2.text(0.5, 0.5, 'No data', ha='center', va='center',
             transform=ax2.transAxes, fontsize=12)
    ax2.set_title('Sex Concordance — No Data', fontsize=12, fontweight='bold')

plt.suptitle('Sex Check QC Report', fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(args.plot, dpi=150, bbox_inches='tight')
plt.close()

# ── Print summary ─────────────────────────────────────────────────────────────
print(f"\n{'='*52}")
print(f"Sex Concordance Summary")
print(f"{'='*52}")
print(f"Total samples             : {len(merged)}")
print(f"Samples with valid F-stat : {merged['F'].notna().sum()}")
if has_collected:
    n_c = (merged['CONCORDANCE'] == 'Concordant').sum()
    n_d = (merged['CONCORDANCE'] == 'Discordant').sum()
    n_a = (merged['CONCORDANCE'] == 'Ambiguous/Unknown').sum()
    print(f"Concordant                : {n_c} ({100*n_c/len(merged):.1f}%)")
    print(f"Discordant                : {n_d} ({100*n_d/len(merged):.1f}%)")
    print(f"Ambiguous/Unknown F-stat  : {n_a} ({100*n_a/len(merged):.1f}%)")
print(f"{'='*52}")

if merged['F'].isna().all():
    print("\nNOTE: All F-statistics are NaN.")
    print("This typically means the chrX dataset was monomorphic after QC,")
    print("or no polymorphic variants survived the sex-stratified filters.")