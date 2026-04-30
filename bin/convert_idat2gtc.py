"""
convert_idat2gtc.py
===================
Converts a single sample's idat files to a GTC file using bcftools +idat2gtc.

bcftools +idat2gtc --idats <dir> writes output GTC files into an output
directory named automatically as <barcode>_<position>.gtc. This script
runs the conversion then renames the output to the desired sample name.

Usage:
    python convert_idat2gtc.py \
        --bpm    /path/to/manifest.bpm \
        --egt    /path/to/clusters.egt \
        --idats  /path/to/per_sample_idat_dir \
        --output sample_id.gtc
"""

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--bpm",    required=True)
parser.add_argument("--egt",    required=True)
parser.add_argument("--idats",  required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

bpm    = Path(args.bpm).resolve()
egt    = Path(args.egt).resolve()
idats  = Path(args.idats).resolve()
output = Path(args.output).resolve()   # desired final path e.g. SAMPLE_ID.gtc

for f, label in [(bpm, "BPM"), (egt, "EGT"), (idats, "idat directory")]:
    if not f.exists():
        sys.exit(f"[convert_idat2gtc] ERROR: {label} not found: {f}")

# bcftools +idat2gtc writes GTC files into a directory
# Use a temp output directory alongside the desired output file
gtc_outdir = output.parent / f"_gtc_tmp_{output.stem}"
gtc_outdir.mkdir(parents=True, exist_ok=True)

cmd = [
    "bcftools", "+idat2gtc",
    "--bpm",    str(bpm),
    "--egt",    str(egt),
    "--idats",  str(idats),
    "--output", str(gtc_outdir),
]

print(f"[convert_idat2gtc] idat dir : {idats}", flush=True)
print(f"[convert_idat2gtc] gtc outdir: {gtc_outdir}", flush=True)
print(f"[convert_idat2gtc] Running  : {' '.join(cmd)}\n", flush=True)

result = subprocess.run(cmd, capture_output=True, text=True)

if result.stdout:
    print(result.stdout, flush=True)
if result.stderr:
    print(result.stderr, file=sys.stderr, flush=True)

if result.returncode != 0:
    # Check if the failure was due to a missing idat file
    stderr_text = result.stderr if result.stderr else ""
    if ("No such file or directory" in stderr_text
            or "Could not open" in stderr_text
            or "Error while running linsolve" in stderr_text
            or "linsolve" in stderr_text.lower()):
        print(
            f"[convert_idat2gtc] WARNING: Some samples could not be processed in {idats}. "
            f"Reason may be missing idat files or numerical failure during normalization. "
            f"Skipping affected samples and continuing with available pairs.",
            file=sys.stderr, flush=True
        )
        # Check if any GTCs were produced before the failure
        gtc_files = list(gtc_outdir.glob("*.gtc")) if gtc_outdir.exists() else []
        if not gtc_files:
            import shutil
            shutil.rmtree(str(gtc_outdir), ignore_errors=True)
            sys.exit(
                f"[convert_idat2gtc] ERROR: No GTC files produced and idat files missing. "
                f"This sample group cannot be processed."
            )
        print(f"[convert_idat2gtc] Partial success: {len(gtc_files)} GTC(s) produced.", flush=True)
    else:
        sys.exit(
            f"[convert_idat2gtc] ERROR: bcftools +idat2gtc failed "
            f"(exit code {result.returncode})"
        )

# Find the GTC file written into the temp output directory
gtc_files = list(gtc_outdir.glob("*.gtc"))
if not gtc_files:
    sys.exit(
        f"[convert_idat2gtc] ERROR: No GTC file produced in {gtc_outdir}\n"
        f"  Contents: {list(gtc_outdir.iterdir())}"
    )

if len(gtc_files) > 1:
    print(
        f"[convert_idat2gtc] WARNING: Multiple GTC files found: {gtc_files}\n"
        f"  Using first: {gtc_files[0]}",
        file=sys.stderr
    )

produced_gtc = gtc_files[0]
print(f"[convert_idat2gtc] Produced  : {produced_gtc.name}", flush=True)

# Move/rename to the desired output path
shutil.move(str(produced_gtc), str(output))

# Clean up temp directory
shutil.rmtree(str(gtc_outdir), ignore_errors=True)

if not output.exists():
    sys.exit(f"[convert_idat2gtc] ERROR: Output GTC not found after rename: {output}")

print(f"[convert_idat2gtc] Done → {output}", flush=True)
