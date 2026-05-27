#!/usr/bin/env python3
"""
Memory-efficient pairwise genotype similarity.
Never builds a wide merged matrix. Instead reads two sample TSVs at a time,
joins them on rsid, and computes similarity — keeping peak memory low.
"""
import sys
import itertools
import pandas as pd
from pathlib import Path


def load_sample(tsv_path):
    """Read a single-sample TSV (rsid index + one genotype column)."""
    df = pd.read_csv(tsv_path, sep="\t", index_col="rsid", dtype=str)
    return df


def compare_pair(tsv_a, tsv_b):
    """
    Load two sample TSVs, join on rsid, count matching genotype calls.
    Returns (sample_a_name, sample_b_name, similarity_pct, matches, total).
    """
    a = load_sample(tsv_a)
    b = load_sample(tsv_b)

    sample_a = a.columns[0]
    sample_b = b.columns[0]

    # Inner join: only SNPs called in BOTH samples
    joined = a.join(b, how="inner")

    col_a = joined.iloc[:, 0]
    col_b = joined.iloc[:, 1]

    # Drop rows where either call is NaN
    valid = col_a.notna() & col_b.notna()
    col_a = col_a[valid]
    col_b = col_b[valid]

    total = len(col_a)
    if total == 0:
        return sample_a, sample_b, 0.0, 0, 0

    matches = int((col_a == col_b).sum())
    similarity = (matches / total) * 100
    return sample_a, sample_b, similarity, matches, total


def main():
    if len(sys.argv) < 3:
        print("Usage: compare_files.py <output_file> <sample1.tsv> <sample2.tsv> ...")
        sys.exit(1)

    outfile = sys.argv[1]
    tsv_files = sys.argv[2:]

    print(f"Samples to compare: {len(tsv_files)}")
    for f in tsv_files:
        print(f"  {f}")

    pairs = list(itertools.combinations(tsv_files, 2))
    print(f"\nRunning {len(pairs)} pairwise comparisons...\n")

    results = []
    for i, (tsv_a, tsv_b) in enumerate(pairs, 1):
        s1, s2, similarity, matches, total = compare_pair(tsv_a, tsv_b)
        line = f"{s1} vs {s2}: {similarity:.2f}% ({matches}/{total} SNPs in common called)"
        results.append(line)
        print(f"[{i}/{len(pairs)}] {line}")

    with open(outfile, "w") as f:
        f.write("# Pairwise genotype similarity (aligned on rsID — cross-array safe)\n")
        f.write("# Only SNPs with a call in BOTH samples are counted\n\n")
        for line in results:
            f.write(line + "\n")

    print(f"\nDone. Results written to {outfile}")


if __name__ == "__main__":
    main()
