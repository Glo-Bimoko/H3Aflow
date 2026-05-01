#!/usr/bin/env python
"""
flag_ibd_duplicates.py
======================
Identifies duplicate/identical samples from IBD/KING output.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--genome", help="PLINK1.9 .genome file")
parser.add_argument("--king", help="PLINK2 .king.kin0 file")
parser.add_argument("--pi_hat", type=float, default=0.1875, help="PI_HAT cutoff for relatedness")
parser.add_argument("--out", required=True, help="Output file for all related pairs")
parser.add_argument("--identical", required=True, help="Output file for genetically identical samples")
parser.add_argument("--plot", required=True, help="Output plot file")
args = parser.parse_args()

# Read input file
if args.genome:
    # PLINK1.9 format - fix the regex warning
    df = pd.read_csv(args.genome, sep=r'\s+')
    # PI_HAT > cutoff indicates duplicates/first-degree relatives
    df['RELATED'] = df['PI_HAT'] > args.pi_hat
    kinship_col = 'PI_HAT'
elif args.king:
    # PLINK2 KING format
    df = pd.read_csv(args.king, sep=r'\s+')
    # Kinship coefficient > 0.0884 indicates duplicates/1st degree
    df['RELATED'] = df['Kinship'] > 0.0884
    kinship_col = 'Kinship'
else:
    raise ValueError("Either --genome or --king must be provided")

# Find genetically identical samples (PI_HAT ~ 0.5 or Kinship ~ 0.25)
if args.genome:
    identical = df[df['PI_HAT'] > 0.45]  # Monozygotic twins/duplicates
    identical_cutoff = 0.45
else:
    identical = df[df['Kinship'] > 0.177]  # Duplicates in KING
    identical_cutoff = 0.177

# Save results
df.to_csv(args.out, sep="\t", index=False)
identical.to_csv(args.identical, sep="\t", index=False)

# Create plot
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: Distribution of kinship/PI_HAT
ax1 = axes[0]
if args.genome:
    ax1.hist(df['PI_HAT'], bins=50, alpha=0.7, color='steelblue', edgecolor='black')
    ax1.axvline(x=args.pi_hat, color='red', linestyle='--', linewidth=2, 
                label=f'Relatedness cutoff (PI_HAT={args.pi_hat})')
    ax1.axvline(x=identical_cutoff, color='orange', linestyle='--', linewidth=2, 
                label=f'Identical cutoff (PI_HAT>{identical_cutoff})')
    ax1.set_xlabel('PI_HAT (Proportion IBD)', fontsize=12)
else:
    ax1.hist(df['Kinship'], bins=50, alpha=0.7, color='steelblue', edgecolor='black')
    ax1.axvline(x=0.0884, color='red', linestyle='--', linewidth=2, 
                label=f'Relatedness cutoff (Kinship>0.0884)')
    ax1.axvline(x=identical_cutoff, color='orange', linestyle='--', linewidth=2, 
                label=f'Identical cutoff (Kinship>{identical_cutoff})')
    ax1.set_xlabel('Kinship Coefficient', fontsize=12)
ax1.set_ylabel('Count', fontsize=12)
ax1.set_title('Pairwise Relatedness Distribution', fontsize=14, fontweight='bold')
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

# Plot 2: Summary statistics
ax2 = axes[1]
categories = ['Total Pairs', 'Related Pairs', 'Identical Pairs']
values = [len(df), df['RELATED'].sum(), len(identical)]
colors = ['lightblue', 'lightcoral', 'gold']
bars = ax2.bar(categories, values, color=colors, alpha=0.7, edgecolor='black')
ax2.set_ylabel('Count', fontsize=12)
ax2.set_title('IBD Summary Statistics', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3, axis='y')

# Add value labels on bars
for bar, value in zip(bars, values):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height + max(values)*0.01,
             f'{value:,}', ha='center', va='bottom', fontsize=10, fontweight='bold')

plt.tight_layout()
plt.savefig(args.plot, dpi=150, bbox_inches='tight')
plt.close()

# Print summary to console
print(f"\n{'='*50}")
print(f"IBD Analysis Summary")
print(f"{'='*50}")
print(f"Total pairs analyzed: {len(df):,}")
print(f"Related pairs (PI_HAT > {args.pi_hat}): {df['RELATED'].sum():,}")
print(f"Genetically identical pairs (PI_HAT > {identical_cutoff}): {len(identical):,}")
print(f"{'='*50}")

if len(identical) > 0:
    print(f"\nIdentical samples (potential duplicates):")
    print(f"{'-'*50}")
    for idx, row in identical.head(20).iterrows():  # Show first 20
        if args.genome:
            # Convert to string to handle integer IDs
            sample1 = str(row['IID1'])
            sample2 = str(row['IID2'])
            print(f"  {sample1} <-> {sample2}: PI_HAT={row['PI_HAT']:.4f}")
        else:
            sample1 = str(row['IID1'])
            sample2 = str(row['IID2'])
            print(f"  {sample1} <-> {sample2}: Kinship={row['Kinship']:.4f}")
    
    if len(identical) > 20:
        print(f"  ... and {len(identical)-20} more pairs")
    print(f"{'-'*50}")
    
    # Create a list of unique samples that are duplicates
    duplicate_samples = set()
    for _, row in identical.iterrows():
        duplicate_samples.add(str(row['IID1']))  # Convert to string
        duplicate_samples.add(str(row['IID2']))  # Convert to string
    
    print(f"\nSamples involved in duplicate pairs: {len(duplicate_samples)}")
    samples_list = list(duplicate_samples)[:10]
    print(f"  {', '.join(samples_list)}")
    if len(duplicate_samples) > 10:
        print(f"  ... and {len(duplicate_samples)-10} more")
else:
    print(f"\nNo genetically identical pairs found.")