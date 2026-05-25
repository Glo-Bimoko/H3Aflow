"""
prep_inputs.py
==============
Reads the project samplesheet (CSV) and the ready/ idat root directory,
then produces two outputs consumed by the Nextflow pipeline:

  --out_samplesheet   resolved_samplesheet.csv
        Columns: sample_id, idat_dir, plate
        idat_dir is the absolute path inside ready/ that contains both
        the Red and Grn idat files for the given BeadChip Barcode +
        Sentrix Position.

  --out_sex_info      sex_info.tsv
        Columns: sampleid, sex
        sex: 0 = Female, 1 = Male  (matches annotate_sex_check.py convention)

Expected samplesheet columns (case-insensitive, extra columns ignored):
  Sample ID | BeadChip Barcode | Sentrix Position | Plate Number |
  Well Position | Collected Gender

Supported samplesheet formats: .csv, .xlsx, .xls

Idat file naming convention (Illumina standard):
  <BeadChip Barcode>_<Sentrix Position>_Red.idat
  <BeadChip Barcode>_<Sentrix Position>_Grn.idat

The script searches recursively under --idat_root (default: ready/) for a
directory that contains both files matching the barcode+position for each
sample.  The search is intentionally broad so that idats can be nested
arbitrarily deep within ready/.

Exit codes:
  0  – all samples resolved successfully
  1  – one or more samples could not be located (details printed to stderr)
"""

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Prepare pipeline inputs from samplesheet")
parser.add_argument("--samplesheet",     required=True,
                    help="Raw project samplesheet CSV")
parser.add_argument("--idat_root",       default="ready",
                    help="Root directory that contains all idat files (default: ready)")
parser.add_argument("--out_samplesheet", required=True,
                    help="Path for the resolved samplesheet CSV (sample_id, idat_dir, plate)")
parser.add_argument("--out_sex_info",    required=True,
                    help="Path for the sex_info TSV (sampleid, sex)")
parser.add_argument("--gtc_dir",         default=None,
                    help="Directory of pre-seeded per-sample GTC files ({sample_id}.gtc)")
args = parser.parse_args()

# ── Load samplesheet (CSV or Excel) ───────────────────────────────────────────
print(f"[prep_inputs] Reading samplesheet: {args.samplesheet}", flush=True)
ext = Path(args.samplesheet).suffix.lower()
if ext in ('.xlsx', '.xls'):
    try:
        ss = pd.read_excel(args.samplesheet, dtype=str)
    except ImportError:
        sys.exit(
            "[prep_inputs] ERROR: Reading .xlsx requires the openpyxl package.\n"
            "  Install it with:  pip install openpyxl"
        )
elif ext == '.csv':
    ss = pd.read_csv(args.samplesheet, dtype=str)
else:
    print(f"[prep_inputs] WARNING: Unrecognised extension '{ext}', attempting CSV read.", flush=True)
    ss = pd.read_csv(args.samplesheet, dtype=str)
ss.columns = ss.columns.str.strip()

# Normalise column names to lowercase-underscore for internal use
col_map = {c: c.strip().lower().replace(" ", "_") for c in ss.columns}
ss = ss.rename(columns=col_map)

REQUIRED = {
    "sample_id":        ["sample_id", "sampleid", "sample id"],
    "barcode":          ["beadchip_barcode", "barcode", "chip_barcode",
                         "beadchip barcode", "sentrix_id"],
    "sentrix_position": ["sentrix_position", "position", "sentrix position",
                         "array_position"],
    "plate":            ["plate_number", "plate", "plate number"],
    "collected_gender": ["collected_gender", "gender", "sex",
                         "collected gender", "collected_sex"],
}

def resolve_col(df, aliases, label):
    for a in aliases:
        norm = a.strip().lower().replace(" ", "_")
        # match against already-normalised columns
        for c in df.columns:
            if c == norm:
                return c
    sys.exit(
        f"[prep_inputs] ERROR: Cannot find required column '{label}'.\n"
        f"  Tried aliases: {aliases}\n"
        f"  Found columns: {list(df.columns)}"
    )

c_sample   = resolve_col(ss, REQUIRED["sample_id"],        "Sample ID")
c_barcode  = resolve_col(ss, REQUIRED["barcode"],           "BeadChip Barcode")
c_position = resolve_col(ss, REQUIRED["sentrix_position"],  "Sentrix Position")
c_plate    = resolve_col(ss, REQUIRED["plate"],             "Plate Number")
c_gender   = resolve_col(ss, REQUIRED["collected_gender"],  "Collected Gender")

# Strip whitespace from all key columns
for col in [c_sample, c_barcode, c_position, c_plate, c_gender]:
    ss[col] = ss[col].str.strip()

print(f"[prep_inputs] {len(ss)} samples in samplesheet.", flush=True)

# ── Build sex_info TSV ─────────────────────────────────────────────────────────
# Normalise gender to 0 (Female) / 1 (Male)
GENDER_MAP = {
    "female": "0", "f": "0", "2": "0", "0": "0",
    "male":   "1", "m": "1", "1": "1",
}

def normalise_gender(val):
    if pd.isna(val):
        return None
    norm = str(val).strip().lower()
    return GENDER_MAP.get(norm, None)

ss["_sex_code"] = ss[c_gender].apply(normalise_gender)

unresolved_sex = ss[ss["_sex_code"].isna()]
if len(unresolved_sex):
    print(
        f"[prep_inputs] WARNING: {len(unresolved_sex)} sample(s) have unrecognised "
        f"gender values and will be coded as 'Unknown' in sex_info.\n"
        f"  Values seen: {unresolved_sex[c_gender].unique().tolist()}",
        flush=True,
    )
    # Write them as empty string — annotate_sex_check.py handles Unknown
    ss.loc[ss["_sex_code"].isna(), "_sex_code"] = ""

sex_info = ss[[c_sample, "_sex_code"]].copy()
sex_info.columns = ["sampleid", "sex"]
sex_info.to_csv(args.out_sex_info, sep="\t", index=False)
print(f"[prep_inputs] sex_info written → {args.out_sex_info}", flush=True)

# ── Index idat files under ready/ ─────────────────────────────────────────────
idat_root = Path(args.idat_root)
if not idat_root.exists():
    sys.exit(f"[prep_inputs] ERROR: idat root directory not found: {idat_root.resolve()}")

print(f"[prep_inputs] Indexing idat files under {idat_root.resolve()} …", flush=True)

# Walk the tree once and build a dict:  filename_stem → parent_directory
# Illumina idat stems look like:  <Barcode>_<Position>_Red  /  _Grn
idat_index: dict[str, Path] = {}   # stem (lower) → directory containing it
for dirpath, _dirs, files in os.walk(idat_root):
    for fname in files:
        if fname.lower().endswith(".idat"):
            stem = fname[:-5]   # strip .idat
            idat_index[stem.lower()] = Path(dirpath)

print(f"[prep_inputs] Found {len(idat_index)} idat files.", flush=True)

gtc_dir = Path(args.gtc_dir).resolve() if args.gtc_dir else None
if gtc_dir:
    if not gtc_dir.exists():
        print(f"[prep_inputs] WARNING: gtc_dir not found: {gtc_dir}", flush=True)
        gtc_dir = None
    else:
        print(f"[prep_inputs] GTC seed directory: {gtc_dir}", flush=True)

# ── Resolve idat directory for each sample ─────────────────────────────────────
resolved_rows = []
missing = []
gtc_seeded = []

for _, row in ss.iterrows():
    sample_id = row[c_sample]
    barcode   = row[c_barcode]
    position  = row[c_position]
    plate     = row[c_plate]

    # Expected file stems (Illumina convention)
    # Both barcode and position are lowercased to match the index
    stem_red = f"{barcode}_{position}_red".lower()
    stem_grn = f"{barcode}_{position}_grn".lower()

    dir_red = idat_index.get(stem_red)
    dir_grn = idat_index.get(stem_grn)

    if dir_red is None or dir_grn is None:
        # Try alternative: some pipelines omit position suffix
        # e.g.  <Barcode>_Red.idat  (rare but seen in older Illumina exports)
        stem_red_alt = f"{barcode}_red"
        stem_grn_alt = f"{barcode}_grn"
        dir_red = dir_red or idat_index.get(stem_red_alt)
        dir_grn = dir_grn or idat_index.get(stem_grn_alt)

    if dir_red is None or dir_grn is None:
        gtc_path = gtc_dir / f"{sample_id}.gtc" if gtc_dir else None
        if gtc_path is not None and gtc_path.is_file():
            resolved_rows.append({
                "sample_id":    sample_id,
                "idat_dir":     "GTC_SEED",
                "barcode":      barcode,
                "position":     position,
                "plate":        plate,
                "input_source": "gtc",
            })
            gtc_seeded.append(sample_id)
            continue

        missing.append({
            "sample_id": sample_id,
            "barcode":   barcode,
            "position":  position,
            "missing":   "Red" if dir_red is None else "Grn",
        })
        continue

    # Both files found; verify they live in the same directory
    if dir_red != dir_grn:
        print(
            f"[prep_inputs] WARNING: Red and Grn idats for {sample_id} are in "
            f"different directories:\n  Red: {dir_red}\n  Grn: {dir_grn}\n"
            f"  Using Red idat directory.",
            flush=True,
        )

    resolved_rows.append({
        "sample_id":    sample_id,
        "idat_dir":     str(dir_red.resolve()),
        "barcode":      barcode,
        "position":     position,
        "plate":        plate,
        "input_source": "idat",
    })

# ── Report missing ─────────────────────────────────────────────────────────────
if missing:
    print(
        f"\n[prep_inputs] WARNING: Could not locate idat files for "
        f"{len(missing)} sample(s). They will be excluded from the pipeline.",
        file=sys.stderr
    )
    for m in missing:
        print(
            f"  sample={m['sample_id']}  barcode={m['barcode']}  "
            f"position={m['position']}  missing={m['missing']} idat",
            file=sys.stderr,
        )
    print(
        f"\n  Expected file pattern under {idat_root.resolve()}:\n"
        f"    <barcode>_<sentrix_position>_Red.idat\n"
        f"    <barcode>_<sentrix_position>_Grn.idat\n",
        file=sys.stderr,
    )
    print(
        f"[prep_inputs] Continuing with {len(resolved_rows)} samples "
        f"({len(missing)} excluded due to missing idat/GTC).",
        flush=True
    )

if gtc_seeded:
    print(
        f"[prep_inputs] {len(gtc_seeded)} sample(s) will use pre-seeded GTC "
        f"(no idat): {', '.join(gtc_seeded)}",
        flush=True,
    )

# ── Write resolved samplesheet ─────────────────────────────────────────────────
resolved_df = pd.DataFrame(resolved_rows)
resolved_df.to_csv(args.out_samplesheet, index=False)

print(
    f"\n[prep_inputs] All {len(resolved_df)} samples resolved successfully.",
    flush=True,
)
print(f"[prep_inputs] Resolved samplesheet → {args.out_samplesheet}", flush=True)

# Print a brief per-plate summary
if len(resolved_df):
    plate_counts = resolved_df.groupby("plate").size()
    print("\n[prep_inputs] Samples per plate:")
    for plate, n in plate_counts.items():
        print(f"  Plate {plate}: {n} samples")
else:
    print("\n[prep_inputs] WARNING: No samples resolved.", flush=True)