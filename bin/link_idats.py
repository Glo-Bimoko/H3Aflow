"""
link_idats.py
=============
Creates a clean per-sample directory containing symlinks to the Red and Grn
idat files for a given sample. This isolates each sample's idats so that
downstream tools (idat2gtc) receive a directory containing only the files
they need.

Usage:
    python link_idats.py --source /path/to/idat_dir --dest linked_idats/SAMPLE_ID

The source directory is searched for any .idat files and symlinked into dest/.
"""

import argparse
import os
import sys
from pathlib import Path

parser = argparse.ArgumentParser(description="Symlink idat files into a clean directory")
parser.add_argument("--source", required=True,
                    help="Directory containing the idat files for this sample")
parser.add_argument("--dest",   required=True,
                    help="Destination directory to create with symlinks")
args = parser.parse_args()

source = Path(args.source).resolve()
dest   = Path(args.dest).resolve()

if not source.exists():
    sys.exit(
        f"[link_idats] ERROR: Source directory does not exist: {source}\n"
        f"  Check that the idat_dir in resolved_samplesheet.csv is correct."
    )

# Collect all idat files in the source directory
idat_files = list(source.glob("*.idat"))

if not idat_files:
    # Try case-insensitive search (some systems have .IDAT)
    idat_files = [f for f in source.iterdir()
                  if f.suffix.lower() == ".idat"]

if not idat_files:
    sys.exit(
        f"[link_idats] ERROR: No .idat files found in: {source}\n"
        f"  Contents: {[f.name for f in source.iterdir()]}"
    )

# Create destination directory
dest.mkdir(parents=True, exist_ok=True)

# Create symlinks
linked = []
for idat in idat_files:
    link = dest / idat.name
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(idat)
    linked.append(idat.name)

print(f"[link_idats] Linked {len(linked)} idat file(s) from {source} → {dest}",
      flush=True)
for name in sorted(linked):
    print(f"  {name}", flush=True)
