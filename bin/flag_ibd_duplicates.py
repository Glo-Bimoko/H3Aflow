#!/usr/bin/env python3
"""
Simplified flag_ibd_duplicates.py
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--genome", required=True, help="PLINK1.9 .genome file")
parser.add_argument("--pi_hat", type=float, default=0.1875)
parser.add_argument("--out", required=True)
parser.add_argument("--identical", required=True)
parser.add_argument("--plot", required=True)
args = parser.parse_args()

# Read the IBD results
df = pd.read_csv(args.genome, sep=r'\s+')

# Find related and identical pairs
df['RELATED'] = df['PI_HAT'] > args.pi_hat
identical = df[df['PI_HAT'] > 0.45]

# Save outputs
df.to_csv(args.out, sep="\t", index=False)
identical.to_csv(args.identical, sep="\t", index=False)

# Create simple plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Histogram of PI_HAT
ax1.hist(df['PI_HAT'], bins=50, alpha=0.7, color='steelblue', edgecolor='black')
ax1.axvline(x=args.pi_hat, color='red', linestyle='--', label=f'Relatedness cutoff (PI_HAT={args.pi_hat})')
ax1.axvline(x=0.45, color='orange', linestyle='--', label='Identical cutoff (PI_HAT>0.45)')
ax1.set_xlabel('PI_HAT')
ax1.set_ylabel('Count')
ax1.set_title('Pairwise Relatedness Distribution')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Summary bar plot
categories = ['Total Pairs', 'Related Pairs', 'Identical Pairs']
values = [len(df), df['RELATED'].sum(), len(identical)]
colors = ['lightblue', 'lightcoral', 'gold']
ax2.bar(categories, values, color=colors, alpha=0.7, edgecolor='black')
ax2.set_ylabel('Count')
ax2.set_title('IBD Summary')
ax2.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(args.plot, dpi=150, bbox_inches='tight')
plt.close()

# Print summary
print(f"\n{'='*50}")
print(f"IBD Analysis Summary")
print(f"{'='*50}")
print(f"Total pairs analyzed: {len(df):,}")
print(f"Related pairs (PI_HAT > {args.pi_hat}): {df['RELATED'].sum():,}")
print(f"Genetically identical pairs: {len(identical):,}")
if len(identical) > 0:
    print(f"\nFirst 5 identical pairs:")
    for _, row in identical.head(5).iterrows():
        print(f"  {row['IID1']} <-> {row['IID2']}: PI_HAT={row['PI_HAT']:.4f}")
print(f"{'='*50}")