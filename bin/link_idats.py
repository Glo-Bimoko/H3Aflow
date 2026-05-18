#!/usr/bin/env python3
"""
link_idats.py
=============
Creates a clean per-sample directory containing symlinks to the Red and Grn
idat files for a given sample. Filters by barcode+position so that only the
two files belonging to this sample are linked — not the entire plate directory.

Usage:
    python link_idats.py \
        --source   /path/to/plate_idat_dir \
        --dest     linked_idats/SAMPLE_ID \
        --barcode  205695340025 \
        --position R01C01
"""
import argparse
import sys
from pathlib import Path

parser = argparse.ArgumentParser(description="Symlink per-sample idat files into a clean directory")
parser.add_argument("--source",   required=True,  help="Plate directory containing all idat files")
parser.add_argument("--dest",     required=True,  help="Destination directory to create with symlinks")
parser.add_argument("--barcode",  required=True,  help="BeadChip barcode for this sample")
parser.add_argument("--position", required=True,  help="Sentrix position for this sample (e.g. R01C01)")
args = parser.parse_args()

source   = Path(args.source).resolve()
dest     = Path(args.dest).resolve()
barcode  = args.barcode.strip()
position = args.position.strip()

if not source.exists():
    sys.exit(
        f"[link_idats] ERROR: Source directory does not exist: {source}\n"
        f"  Check that idat_dir in resolved_samplesheet.csv is correct."
    )

# Build expected filenames for this sample (Illumina convention)
prefix     = f"{barcode}_{position}_"
red_name   = f"{prefix}Red.idat"
grn_name   = f"{prefix}Grn.idat"

# Case-insensitive search across all idats in the plate directory
all_idats  = {f.name.lower(): f for f in source.iterdir() if f.suffix.lower() == ".idat"}

red_file   = all_idats.get(red_name.lower())
grn_file   = all_idats.get(grn_name.lower())

missing = []
if red_file is None:
    missing.append(red_name)
if grn_file is None:
    missing.append(grn_name)

if missing:
    sys.exit(
        f"[link_idats] ERROR: Could not find idat file(s) for "
        f"barcode={barcode} position={position}:\n"
        f"  Missing: {missing}\n"
        f"  Searched in: {source}\n"
        f"  Available files (first 10): {sorted(all_idats.keys())[:10]}"
    )

# Create destination and symlink only these two files
dest.mkdir(parents=True, exist_ok=True)

for src_file in (red_file, grn_file):
    link = dest / src_file.name
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(src_file)

print(
    f"[link_idats] Linked 2 idat file(s) for {barcode}_{position} → {dest}",
    flush=True,
)
print(f"  {red_file.name}", flush=True)
print(f"  {grn_file.name}", flush=True)