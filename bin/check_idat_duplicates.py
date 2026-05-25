#!/usr/bin/env python3
"""
check_idat_duplicates.py
Detect samples that share the same BeadChip barcode + Sentrix position (same
idat source).  Genotype concordance only compares samples in the merged cohort;
this script catches intentional or accidental idat duplicates from the
samplesheet even when only one alias reached QC.
"""
import argparse
import hashlib
import sys
from pathlib import Path

import pandas as pd

PLATE_CANDIDATES = [
    "beadchip_barcode", "barcode", "chip_barcode", "sentrix_barcode",
]
POSITION_CANDIDATES = [
    "sentrix_position", "position", "sentrix_pos",
]
SAMPLE_CANDIDATES = ["sample_id", "sampleid", "sample id"]


def parse_args():
    parser = argparse.ArgumentParser(description="Flag samples sharing the same idat source")
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--fam", required=True, help="Cohort .fam file (post-QC sample list)")
    parser.add_argument("--concordance", default=None,
                        help="pairwise_concordance.tsv (optional, for pairs both in cohort)")
    parser.add_argument("--out", required=True, help="Output TSV report")
    parser.add_argument("--log", default=None, help="Append summary to this log file")
    parser.add_argument("--gtc_dir", default=None,
                        help="Optional GTC directory to detect identical GTC content")
    return parser.parse_args()


def normalize_columns(samplesheet):
    df = pd.read_csv(samplesheet, dtype=str)
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    return df


def pick_column(columns, candidates):
    for name in candidates:
        if name in columns:
            return name
    return None


def load_cohort_iids(fam_path):
    fam = pd.read_csv(fam_path, sep=r"\s+", header=None,
                      names=["FID", "IID", "FATHER", "MOTHER", "SEX", "PHENO"])
    return set(fam["IID"].astype(str).str.strip())


def load_concordance_lookup(concordance_path):
    if not concordance_path or not Path(concordance_path).exists():
        return {}
    conc = pd.read_csv(concordance_path, sep="\t")
    lookup = {}
    for _, row in conc.iterrows():
        a = str(row["SAMPLE_A"]).strip()
        b = str(row["SAMPLE_B"]).strip()
        key = tuple(sorted((a, b)))
        lookup[key] = row.get("CONCORDANCE_PCT")
    return lookup


def build_idat_key(row, barcode_col, position_col):
    return f"{row[barcode_col].strip()}_{row[position_col].strip()}"


def summarize_group(group, sample_col, cohort_iids, concordance_lookup):
    sample_ids = group[sample_col].astype(str).str.strip().tolist()
    in_cohort = [sid for sid in sample_ids if sid in cohort_iids]
    rows = []
    for i in range(len(sample_ids)):
        for j in range(i + 1, len(sample_ids)):
            a, b = sample_ids[i], sample_ids[j]
            key = tuple(sorted((a, b)))
            rows.append({
                "SAMPLE_A": a,
                "SAMPLE_B": b,
                "SAMPLE_A_IN_COHORT": a in cohort_iids,
                "SAMPLE_B_IN_COHORT": b in cohort_iids,
                "BOTH_IN_COHORT": a in cohort_iids and b in cohort_iids,
                "CONCORDANCE_PCT": concordance_lookup.get(key),
                "STATUS": classify_pair(a in cohort_iids, b in cohort_iids, concordance_lookup.get(key)),
            })
    return rows, in_cohort, sample_ids


def classify_pair(a_in, b_in, concordance_pct):
    if a_in and b_in:
        if concordance_pct is None or pd.isna(concordance_pct):
            return "BOTH_IN_COHORT_NO_CONCORDANCE"
        if concordance_pct >= 99.0:
            return "GENOTYPE_DUPLICATE"
        return "BOTH_IN_COHORT_LOW_CONCORDANCE"
    if a_in or b_in:
        return "SINGLE_ALIAS_IN_COHORT"
    return "NEITHER_IN_COHORT"


def append_log(log_path, report, idat_key, sample_ids, in_cohort):
    with open(log_path, "a", encoding="utf-8") as handle:
        handle.write(f"\nIdat source: {idat_key}\n")
        handle.write(f"  Sample IDs in samplesheet : {', '.join(sample_ids)}\n")
        handle.write(f"  Sample IDs in QC cohort   : {', '.join(in_cohort) if in_cohort else '(none)'}\n")
        for _, row in report.iterrows():
            handle.write(
                f"  {row['SAMPLE_A']} vs {row['SAMPLE_B']}: "
                f"{row['STATUS']}"
            )
            if pd.notna(row.get("CONCORDANCE_PCT")):
                handle.write(f" ({row['CONCORDANCE_PCT']:.4f}%)")
            handle.write("\n")


def find_identical_gtc_groups(gtc_dir):
    gtc_path = Path(gtc_dir)
    if not gtc_path.is_dir():
        return []

    by_hash = {}
    for path in sorted(gtc_path.glob("*.gtc")):
        digest = hashlib.md5(path.read_bytes()).hexdigest()
        by_hash.setdefault(digest, []).append(path.stem)

    return [names for names in by_hash.values() if len(names) > 1]


def report_identical_gtc_groups(groups, cohort_iids, concordance_lookup, log_path):
    if not groups:
        return

    with open(log_path, "a", encoding="utf-8") as handle:
        handle.write("\nIdentical GTC file content (byte-for-byte duplicates):\n")
        for names in groups:
            handle.write(f"  {', '.join(names)}\n")
            in_cohort = [name for name in names if name in cohort_iids]
            handle.write(f"    In QC cohort: {', '.join(in_cohort) if in_cohort else '(none)'}\n")
            if len(in_cohort) < 2:
                handle.write(
                    "    NOTE: GTC duplicates only appear in concordance when >=2 "
                    "Sample IDs are in the QC cohort.\n"
                )
                continue
            for i in range(len(in_cohort)):
                for j in range(i + 1, len(in_cohort)):
                    key = tuple(sorted((in_cohort[i], in_cohort[j])))
                    pct = concordance_lookup.get(key)
                    if pct is not None:
                        handle.write(
                            f"    {in_cohort[i]} vs {in_cohort[j]}: "
                            f"genotype concordance {pct:.4f}%\n"
                        )


def main():
    args = parse_args()
    samplesheet = normalize_columns(args.samplesheet)

    sample_col = pick_column(samplesheet.columns, SAMPLE_CANDIDATES)
    barcode_col = pick_column(samplesheet.columns, PLATE_CANDIDATES)
    position_col = pick_column(samplesheet.columns, POSITION_CANDIDATES)
    if not all([sample_col, barcode_col, position_col]):
        sys.exit("[check_idat_duplicates] ERROR: samplesheet missing required columns")

    cohort_iids = load_cohort_iids(args.fam)
    concordance_lookup = load_concordance_lookup(args.concordance)

    samplesheet = samplesheet.copy()
    samplesheet["_idat_key"] = samplesheet.apply(
        lambda row: build_idat_key(row, barcode_col, position_col), axis=1,
    )

    dup_df = samplesheet.groupby("_idat_key").filter(lambda grp: len(grp) > 1)
    grouped = list(dup_df.groupby("_idat_key"))

    all_rows = []
    n_groups = len(grouped)
    for idat_key, group in grouped:
        pair_rows, in_cohort, sample_ids = summarize_group(
            group, sample_col, cohort_iids, concordance_lookup,
        )
        for row in pair_rows:
            row["IDAT_KEY"] = idat_key
            all_rows.append(row)

    report = pd.DataFrame(all_rows)
    if report.empty:
        report = pd.DataFrame(columns=[
            "IDAT_KEY", "SAMPLE_A", "SAMPLE_B", "SAMPLE_A_IN_COHORT",
            "SAMPLE_B_IN_COHORT", "BOTH_IN_COHORT", "CONCORDANCE_PCT", "STATUS",
        ])
    else:
        report = report.sort_values(["IDAT_KEY", "SAMPLE_A", "SAMPLE_B"])

    report.to_csv(args.out, sep="\t", index=False)
    print(f"[check_idat_duplicates] Shared-idat groups: {n_groups}", flush=True)
    print(f"[check_idat_duplicates] Report written → {args.out}", flush=True)

    if args.log and args.gtc_dir:
        gtc_groups = find_identical_gtc_groups(args.gtc_dir)
        if gtc_groups:
            with open(args.log, "a", encoding="utf-8") as handle:
                handle.write("\n" + "=" * 60 + "\n")
                handle.write("Identical GTC content (GTC-level duplicate detection)\n")
                handle.write("=" * 60 + "\n")
            report_identical_gtc_groups(
                gtc_groups, cohort_iids, concordance_lookup, args.log,
            )

    if args.log and n_groups > 0:
        with open(args.log, "a", encoding="utf-8") as handle:
            handle.write("\n" + "=" * 60 + "\n")
            handle.write("Shared idat source (samplesheet duplicate detection)\n")
            handle.write("=" * 60 + "\n")
            handle.write(
                "Genotype concordance only compares samples in the QC cohort.\n"
                "Pairs with STATUS=SINGLE_ALIAS_IN_COHORT share idat files but\n"
                "only one Sample ID passed filtering — concordance cannot compare them.\n\n"
            )
        for idat_key, group in grouped:
            sub = report[report["IDAT_KEY"] == idat_key]
            sample_ids = group[sample_col].astype(str).str.strip().tolist()
            in_cohort = [sid for sid in sample_ids if sid in cohort_iids]
            append_log(args.log, sub, idat_key, sample_ids, in_cohort)

        single_alias = report[report["STATUS"] == "SINGLE_ALIAS_IN_COHORT"]
        if len(single_alias) > 0:
            with open(args.log, "a", encoding="utf-8") as handle:
                handle.write(
                    f"\nWARNING: {len(single_alias)} shared-idat pair(s) have only one "
                    "alias in the QC cohort — include both Sample IDs to validate via concordance.\n"
                )

        with open(args.log, "a", encoding="utf-8") as handle:
            handle.write("=" * 60 + "\n")


if __name__ == "__main__":
    main()
