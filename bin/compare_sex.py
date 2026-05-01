#!/usr/bin/env python
"""
compare_sex.py
==============
Compares PLINK-inferred sex (from chrX F-statistic) with collected sex.
Identifies discordant samples for QC flagging.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument("--plink_sex", required=True, help="PLINK .sexcheck file")
parser.add_argument("--collected_sex", required=True, help="Collected sex info TSV")
parser.add_argument("--out", required=True, help="Output discordance file")
parser.add_argument("--plot", required=True, help="Output plot file")
args = parser.parse_args()

# Read PLINK sex check results
plink_sex = pd.read_csv(args.plink_sex, sep="\s+")
plink_sex['INFERRED_SEX'] = plink_sex['STATUS'].map({
    'OK': plink_sex['F'].apply(lambda x: 'Female' if x < 0.2 else ('Male' if x > 0.8 else 'Unknown')),
    'PROBLEM': 'Unknown'
})

# Read collected sex
collected = pd.read_csv(args.collected_sex, sep="\t")
collected.columns = [c.lower() for c in collected.columns]
collected = collected.rename(columns={'sampleid': 'IID', 'sex': 'COLLECTED_SEX_NUM'})
collected['COLLECTED_SEX'] = collected['COLLECTED_SEX_NUM'].map({0: 'Female', 1: 'Male'})

# Merge and compare
merged = plink_sex.merge(collected[['IID', 'COLLECTED_SEX']], on='IID', how='left')
merged['CONCORDANT'] = merged['INFERRED_SEX'] == merged['COLLECTED_SEX']
merged['DISCORDANCE_TYPE'] = np.where(
    merged['CONCORDANT'], 'Concordant',
    merged['INFERRED_SEX'] + ' vs ' + merged['COLLECTED_SEX']
)

# Save discordant samples
discordant = merged[~merged['CONCORDANT']]
discordant.to_csv(args.out, sep="\t", index=False)

# Create plot
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Plot 1: F-statistic distribution by collected sex
ax1 = axes[0]
for sex, color in [('Female', 'pink'), ('Male', 'lightblue')]:
    subset = merged[merged['COLLECTED_SEX'] == sex]
    ax1.hist(subset['F'], bins=30, alpha=0.5, label=sex, color=color, edgecolor='black')
ax1.axvline(x=0.2, color='red', linestyle='--', label='Female cutoff (0.2)')
ax1.axvline(x=0.8, color='red', linestyle='--', label='Male cutoff (0.8)')
ax1.set_xlabel('F-statistic (chrX inbreeding coefficient)')
ax1.set_ylabel('Count')
ax1.set_title('Sex Inference by F-statistic')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: Concordance summary
ax2 = axes[1]
concordance_counts = merged['DISCORDANCE_TYPE'].value_counts()
colors = ['green' if x == 'Concordant' else 'red' for x in concordance_counts.index]
bars = ax2.bar(range(len(concordance_counts)), concordance_counts.values, color=colors, alpha=0.7)
ax2.set_xticks(range(len(concordance_counts)))
ax2.set_xticklabels(concordance_counts.index, rotation=45, ha='right')
ax2.set_ylabel('Count')
ax2.set_title('Sex Concordance')
ax2.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(args.plot, dpi=150, bbox_inches='tight')
plt.close()

# Print summary
print(f"\nSex Concordance Summary:")
print(f"  Total samples: {len(merged)}")
print(f"  Concordant: {merged['CONCORDANT'].sum()} ({merged['CONCORDANT'].sum()/len(merged)*100:.1f}%)")
print(f"  Discordant: {(~merged['CONCORDANT']).sum()} ({(~merged['CONCORDANT']).sum()/len(merged)*100:.1f}%)")
if len(discordant) > 0:
    print(f"\nDiscordant samples ({len(discordant)}):")
    for _, row in discordant.iterrows():
        print(f"  {row['IID']}: Inferred={row['INFERRED_SEX']}, Collected={row['COLLECTED_SEX']} (F={row['F']:.3f})")