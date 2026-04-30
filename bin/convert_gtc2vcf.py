"""
convert_gtc2vcf.py
==================
Converts a set of GTC files for one plate to a normalised, sorted BCF using
bcftools +gtc2vcf, then produces a TSV of raw intensities via bcftools query.

Usage:
    python convert_gtc2vcf.py \
        --bpm       /path/to/manifest.bpm          \
        --egt       /path/to/clusters.egt           \
        --gtcs      gtc_list.txt                    \
        --fasta     /path/to/reference.fa           \
        --outprefix Plate_1

Outputs:
    <outprefix>.bcf       normalised, sorted BCF (all samples for the plate)
    <outprefix>.bcf.csi   BCF index
    <outprefix>.tsv       raw X/Y intensity table (SAMPLE, CHR, POS, NORMX, NORMY)

Requires:
    bcftools >= 1.11 with gtc2vcf plugin
    BCFTOOLS_PLUGINS env var pointing to the plugin directory
"""

import argparse
import subprocess
import sys
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--bpm",       required=True)
parser.add_argument("--egt",       required=True)
parser.add_argument("--gtcs",      required=True, help="File listing GTC paths, one per line")
parser.add_argument("--fasta",     required=True)
parser.add_argument("--outprefix", required=True)
args = parser.parse_args()

bpm       = Path(args.bpm).resolve()
egt       = Path(args.egt).resolve()
gtcs_list = Path(args.gtcs)
fasta     = Path(args.fasta).resolve()
prefix    = args.outprefix

bcf_out   = Path(f"{prefix}.bcf")
tsv_out   = Path(f"{prefix}.tsv")

# ── Validate inputs ────────────────────────────────────────────────────────────
for f, label in [(bpm, "BPM"), (egt, "EGT"), (fasta, "FASTA"), (gtcs_list, "GTC list")]:
    if not f.exists():
        sys.exit(f"[convert_gtc2vcf] ERROR: {label} not found: {f}")

gtc_files = [l.strip() for l in gtcs_list.read_text().splitlines() if l.strip()]
if not gtc_files:
    sys.exit(f"[convert_gtc2vcf] ERROR: GTC list is empty: {gtcs_list}")

print(f"[convert_gtc2vcf] Plate prefix : {prefix}", flush=True)
print(f"[convert_gtc2vcf] GTC files    : {len(gtc_files)}", flush=True)
print(f"[convert_gtc2vcf] BPM          : {bpm}", flush=True)
print(f"[convert_gtc2vcf] EGT          : {egt}", flush=True)
print(f"[convert_gtc2vcf] FASTA        : {fasta}", flush=True)

# ── Step 1: GTC → unsorted BCF via bcftools +gtc2vcf ─────────────────────────
unsorted_bcf = f"{prefix}_unsorted.bcf"

gtc2vcf_cmd = (
    ["bcftools", "+gtc2vcf",
     "--bpm",    str(bpm),
     "--egt",    str(egt),
     "--fasta",  str(fasta),
     "--output", unsorted_bcf,
     "--output-type", "b",
     "--no-version"]
    + gtc_files
)

print(f"\n[convert_gtc2vcf] Step 1: GTC → BCF", flush=True)
print(f"  {' '.join(gtc2vcf_cmd[:8])} ... [{len(gtc_files)} GTC files]", flush=True)

result = subprocess.run(gtc2vcf_cmd, capture_output=True, text=True)
if result.stdout:
    print(result.stdout, flush=True)
if result.stderr:
    print(result.stderr, file=sys.stderr, flush=True)
if result.returncode != 0:
    sys.exit(f"[convert_gtc2vcf] ERROR: bcftools +gtc2vcf failed (exit {result.returncode})")

# ── Step 2: Sort and index the BCF ────────────────────────────────────────────
print(f"\n[convert_gtc2vcf] Step 2: Sort BCF", flush=True)

sort_cmd = [
    "bcftools", "sort",
    "--output",      str(bcf_out),
    "--output-type", "b",
    "--temp-dir",    ".",
    unsorted_bcf
]

result = subprocess.run(sort_cmd, capture_output=True, text=True)
if result.stderr:
    print(result.stderr, file=sys.stderr, flush=True)
if result.returncode != 0:
    sys.exit(f"[convert_gtc2vcf] ERROR: bcftools sort failed (exit {result.returncode})")

# Remove unsorted intermediate
Path(unsorted_bcf).unlink(missing_ok=True)

# Index
index_cmd = ["bcftools", "index", str(bcf_out)]
result = subprocess.run(index_cmd, capture_output=True, text=True)
if result.returncode != 0:
    sys.exit(f"[convert_gtc2vcf] ERROR: bcftools index failed (exit {result.returncode})")

print(f"[convert_gtc2vcf] BCF written  : {bcf_out}", flush=True)

# ── Step 3: Extract X/Y intensities to TSV ────────────────────────────────────
print(f"\n[convert_gtc2vcf] Step 3: Extract XY intensities → TSV", flush=True)

# Format string extracts per-sample NORMX and NORMY alongside site info
query_cmd = [
    "bcftools", "query",
    "--format", "[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%NORMX\t%NORMY\n]",
    str(bcf_out)
]

with open(tsv_out, "w") as fh:
    fh.write("SAMPLE_ID\tCHR\tPOS\tREF\tALT\tNORMX\tNORMY\n")
    result = subprocess.run(query_cmd, stdout=fh, stderr=subprocess.PIPE, text=True)

if result.stderr:
    print(result.stderr, file=sys.stderr, flush=True)
if result.returncode != 0:
    # TSV is nice-to-have; warn but don't fail the whole process
    print(
        f"[convert_gtc2vcf] WARNING: TSV extraction failed (exit {result.returncode}). "
        f"Downstream XY intensity QC will be skipped for this plate.",
        file=sys.stderr, flush=True
    )
    # Write an empty TSV so Nextflow output pattern is satisfied
    with open(tsv_out, "w") as fh:
        fh.write("SAMPLE_ID\tCHR\tPOS\tREF\tALT\tNORMX\tNORMY\n")
else:
    print(f"[convert_gtc2vcf] TSV written  : {tsv_out}", flush=True)

print(f"\n[convert_gtc2vcf] Done — plate {prefix}", flush=True)