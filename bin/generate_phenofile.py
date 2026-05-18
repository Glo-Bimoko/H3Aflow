#!/usr/bin/env python3
"""
generate_phenofile.py
=====================
Generates a PLINK-format phenotype file (.phe) from the resolved samplesheet
produced by prep_inputs.py.

PLINK .phe format (space-delimited, no header required but we include one):
  FID  IID  PHENO
  FID and IID are both set to sample_id (no family structure assumed).
  PHENO is set to -9 (missing) for all samples at this stage — actual
  phenotype values are added later by the analyst for association testing.

Input (resolved_samplesheet.csv columns):
  sample_id, idat_dir, plate

Output:
  sample.phe  – PLINK phenotype file
"""

import argparse
import sys
import pandas as pd

parser = argparse.ArgumentParser(description="Generate PLINK phenotype file")
parser.add_argument("--samplesheet", required=True,
                    help="Resolved samplesheet CSV (sample_id, idat_dir, plate)")
parser.add_argument("--out", required=True,
                    help="Output .phe file path")
args = parser.parse_args()

print(f"[generate_phenofile] Reading: {args.samplesheet}", flush=True)
ss = pd.read_csv(args.samplesheet, dtype=str)
ss.columns = ss.columns.str.strip().str.lower()

# Resolve sample_id column name tolerantly
ID_ALIASES = ["sample_id", "sampleid", "sample id", "id"]
col_id = None
for a in ID_ALIASES:
    if a in ss.columns:
        col_id = a
        break

if col_id is None:
    sys.exit(
        f"[generate_phenofile] ERROR: Cannot find sample ID column.\n"
        f"  Tried: {ID_ALIASES}\n"
        f"  Found: {list(ss.columns)}"
    )

ss[col_id] = ss[col_id].str.strip()

# Build phenofile: FID = IID = sample_id, PHENO = -9 (missing)
phe = pd.DataFrame({
    "FID":   ss[col_id],
    "IID":   ss[col_id],
    "PHENO": "-9",
})

phe.to_csv(args.out, sep=" ", index=False)

print(f"[generate_phenofile] Written {len(phe)} samples → {args.out}", flush=True)
print(f"[generate_phenofile] NOTE: PHENO column set to -9 (missing) for all samples.")
print(f"[generate_phenofile] Update the PHENO column with real phenotype values")
print(f"[generate_phenofile] before running association testing.")
