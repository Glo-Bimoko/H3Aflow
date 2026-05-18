#!/usr/bin/env python3
"""
compare_sex.py - Handles monomorphic data case
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--plink_sex", required=True)
parser.add_argument("--collected_sex", required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--plot", required=True)
args = parser.parse_args()

# Read PLINK sex check results
plink_sex = pd.read_csv(args.plink_sex, sep=r'\s+')

print(f"Loaded PLINK sex check with {len(plink_sex)} samples")
print(f"Columns: {plink_sex.columns.tolist()}")
print(f"STATUS value counts: {plink_sex['STATUS'].value_counts().to_dict() if 'STATUS' in plink_sex.columns else 'No STATUS column'}")

# Check if F column exists and has valid values
if 'F' not in plink_sex.columns:
    print("Warning: No F-statistic column found in PLINK output")
    plink_sex['F'] = np.nan

# Read collected sex
try:
    collected = pd.read_csv(args.collected_sex, sep="\t")
    collected.columns = [c.lower() for c in collected.columns]
    collected = collected.rename(columns={'sampleid': 'IID', 'sex': 'COLLECTED_SEX_NUM'})
    collected['COLLECTED_SEX'] = collected['COLLECTED_SEX_NUM'].map({0: 'Female', 1: 'Male'}).fillna('Unknown')
    print(f"Loaded collected sex for {len(collected)} samples")
    has_collected = True
except Exception as e:
    print(f"Warning: Could not read collected sex file: {e}")
    has_collected = False
    collected = pd.DataFrame({'IID': plink_sex['IID'], 'COLLECTED_SEX': 'Unknown'})

# Merge data
merged = plink_sex.merge(collected[['IID', 'COLLECTED_SEX']], on='IID', how='left')

# Handle case where all F values are NaN
if merged['F'].isna().all():
    print("All F-statistics are NaN - data is monomorphic (all samples identical)")
    # Assign F values based on collected sex or default
    merged['INFERRED_SEX'] = 'Unknown'
    merged['F'] = np.nan
    merged['STATUS'] = 'PROBLEM'
else:
    # Map F-statistic to inferred sex
    def infer_sex(f_val):
        if pd.isna(f_val):
            return 'Unknown'
        elif f_val < 0.2:
            return 'Female'
        elif f_val > 0.8:
            return 'Male'
        else:
            return 'Unknown'
    
    merged['INFERRED_SEX'] = merged['F'].apply(infer_sex)

# Determine concordance
merged['CONCORDANT'] = merged['INFERRED_SEX'] == merged['COLLECTED_SEX']
merged['DISCORDANCE_TYPE'] = np.where(
    merged['CONCORDANT'], 'Concordant',
    merged['INFERRED_SEX'] + ' vs ' + merged['COLLECTED_SEX']
)

# Save results
discordant = merged[~merged['CONCORDANT']]
discordant.to_csv(args.out, sep="\t", index=False)
merged.to_csv(args.out.replace('.txt', '_all.tsv'), sep="\t", index=False)

# Create plot
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Plot 1: F-statistic distribution
ax1 = axes[0]
if merged['F'].notna().any():
    # Filter out NaN values for plotting
    valid_f = merged['F'].dropna()
    if len(valid_f) > 0:
        ax1.hist(valid_f, bins=30, alpha=0.7, color='steelblue', edgecolor='black')
        ax1.axvline(x=0.2, color='red', linestyle='--', linewidth=2, label='Female cutoff (0.2)')
        ax1.axvline(x=0.8, color='red', linestyle='--', linewidth=2, label='Male cutoff (0.8)')
        ax1.set_xlabel('F-statistic (chrX inbreeding coefficient)', fontsize=11)
        ax1.set_ylabel('Count', fontsize=11)
        ax1.set_title('Sex Inference by F-statistic', fontsize=12, fontweight='bold')
        ax1.legend(fontsize=9)
        ax1.grid(True, alpha=0.3)
    else:
        ax1.text(0.5, 0.5, 'No valid F-statistics to display\nAll samples are genetically identical', 
                 ha='center', va='center', transform=ax1.transAxes, fontsize=12)
        ax1.set_title('Sex Inference - Monomorphic Data', fontsize=12, fontweight='bold')
else:
    ax1.text(0.5, 0.5, 'No F-statistics available\nAll samples appear genetically identical', 
             ha='center', va='center', transform=ax1.transAxes, fontsize=12)
    ax1.set_title('Sex Inference - Monomorphic Data', fontsize=12, fontweight='bold')

# Plot 2: Concordance summary
ax2 = axes[1]
if has_collected and len(merged) > 0:
    concordance_counts = merged['DISCORDANCE_TYPE'].value_counts()
    colors = ['green' if x == 'Concordant' else 'red' for x in concordance_counts.index]
    bars = ax2.bar(range(len(concordance_counts)), concordance_counts.values, color=colors, alpha=0.7, edgecolor='black')
    ax2.set_xticks(range(len(concordance_counts)))
    ax2.set_xticklabels(concordance_counts.index, rotation=45, ha='right', fontsize=9)
    ax2.set_ylabel('Count', fontsize=11)
    ax2.set_title('Sex Concordance', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for bar, value in zip(bars, concordance_counts.values):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                 f'{value}', ha='center', va='bottom', fontsize=10, fontweight='bold')
else:
    ax2.text(0.5, 0.5, 'No collected sex data available', 
             ha='center', va='center', transform=ax2.transAxes, fontsize=12)
    ax2.set_title('Sex Concordance - Missing Data', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig(args.plot, dpi=150, bbox_inches='tight')
plt.close()

# Print summary
print(f"\n{'='*50}")
print(f"Sex Concordance Summary")
print(f"{'='*50}")
print(f"Total samples: {len(merged)}")
print(f"Samples with valid F-statistics: {merged['F'].notna().sum()}")
print(f"Samples with PROBLEM status: {(merged['STATUS'] == 'PROBLEM').sum() if 'STATUS' in merged.columns else len(merged)}")
if has_collected:
    print(f"\nConcordant: {merged['CONCORDANT'].sum()} ({merged['CONCORDANT'].sum()/len(merged)*100:.1f}%)")
    print(f"Discordant: {(~merged['CONCORDANT']).sum()} ({(~merged['CONCORDANT']).sum()/len(merged)*100:.1f}%)")
print(f"{'='*50}")

# Explanation of results
if merged['F'].isna().all():
    print(f"\nNOTE: All samples have NaN F-statistics.")
    print(f"This typically occurs when all samples are genetically identical,")
    print(f"or when there is no variation on the X chromosome in your dataset.")
    print(f"For sex checking to work properly, you need genetic variation.")