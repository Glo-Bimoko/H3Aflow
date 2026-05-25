#!/usr/bin/env python3
"""
pairwise_concordance.py
=======================
Memory-efficient all-vs-all pairwise genotype concordance from a PLINK
.traw file (produced by plink --recode A-transpose).

Design
------
The .traw file has one ROW per SNP and one COLUMN per sample, so we can
stream it in fixed-size chunks of SNPs.  For each chunk we:
  1. Parse into a (chunk_size × N_samples) uint8 array (0/1/2/255=missing).
  2. For every pair (i, j) add to two running accumulators:
       - n_called[i,j]  : SNPs where BOTH samples have a valid call
       - n_match[i,j]   : SNPs where both calls are identical
  3. Discard the chunk and move to the next.

Peak memory = chunk_size × N_samples (the chunk array) +
              N_samples × N_samples × 2 (the accumulators).
With chunk=50,000 SNPs and 148 samples:
  chunk  : 50,000 × 148 × 1 byte ≈  7 MB
  accum  : 148 × 148 × 8 bytes × 2 ≈ 350 KB
  Total  : < 10 MB regardless of total SNP count.

Output
------
1. pairwise_concordance.tsv  — one row per pair, tab-separated:
     SAMPLE_A  SAMPLE_B  N_CALLED  N_MATCH  CONCORDANCE_PCT

2. pairwise_concordance.log  — summary statistics and a flag list of
   pairs whose concordance exceeds --dup-threshold (default 90.0 %).

Usage
-----
  pairwise_concordance.py \\
      --traw    cohort_snpqc.traw \\
      --out     pairwise_concordance.tsv \\
      --chunk   50000 \\
      --dup-threshold 90.0
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--traw",          required=True,  help=".traw file from plink --recode A-transpose")
    p.add_argument("--out",           required=True,  help="Output TSV path")
    p.add_argument("--chunk",         type=int, default=50_000,
                   help="SNPs to process per chunk (default: 50000)")
    p.add_argument("--dup-threshold", type=float, default=90.0,
                   help="Concordance %% above which a pair is flagged as a "
                        "likely duplicate (default: 90.0)")
    return p.parse_args()


# ---------------------------------------------------------------------------
# .traw header parsing
# ---------------------------------------------------------------------------
TRAW_META_COLS = 6          # CHR  SNP  (C)M  POS  COUNTED  ALT
MISSING_CALL   = 255        # sentinel for NA / non-numeric

def read_traw_header(traw_path):
    """
    Return (sample_ids, n_meta_cols).
    .traw header line: CHR SNP (C)M POS COUNTED ALT <IID1> <IID2> ...
    Sample IDs in .traw are encoded as FID_IID when FID != 0, or just IID.
    """
    with open(traw_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    sample_cols = header[TRAW_META_COLS:]
    # Strip FID prefix (FID_IID → IID) — consistent with --double-id convention
    sample_ids = [s.split("_", 1)[-1] if "_" in s else s for s in sample_cols]
    return sample_ids


# ---------------------------------------------------------------------------
# Chunk iterator
# ---------------------------------------------------------------------------
def iter_chunks(traw_path, chunk_size, n_samples):
    """
    Yield numpy arrays of shape (rows_in_chunk, n_samples), dtype uint8.
    Missing / non-numeric calls → MISSING_CALL (255).
    Skips the header row.
    """
    buf = np.empty((chunk_size, n_samples), dtype=np.uint8)
    row_in_chunk = 0

    with open(traw_path) as fh:
        fh.readline()                      # skip header
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            calls  = fields[TRAW_META_COLS:]
            for j, c in enumerate(calls):
                if c == "NA" or c == "":
                    buf[row_in_chunk, j] = MISSING_CALL
                else:
                    try:
                        v = int(c)
                        buf[row_in_chunk, j] = v if 0 <= v <= 2 else MISSING_CALL
                    except ValueError:
                        buf[row_in_chunk, j] = MISSING_CALL

            row_in_chunk += 1
            if row_in_chunk == chunk_size:
                yield buf[:row_in_chunk]
                row_in_chunk = 0

    if row_in_chunk > 0:
        yield buf[:row_in_chunk]


# ---------------------------------------------------------------------------
# Core accumulation
# ---------------------------------------------------------------------------
def accumulate(chunk, n_called, n_match):
    """
    Update n_called[i,j] and n_match[i,j] for all pairs using a chunk
    of shape (S, N) where S=SNPs in chunk, N=samples.

    Vectorised over SNPs: for each SNP row we do an outer-product-style
    update over all N*(N-1)/2 pairs without an explicit Python loop over
    pairs.  This keeps the inner loop in numpy C code.
    """
    # valid[s, i] = True if sample i has a non-missing call at SNP s
    valid = (chunk != MISSING_CALL)                      # (S, N) bool

    # For every pair (i, j) with i < j:
    #   both_called[s] = valid[s,i] & valid[s,j]
    #   matched[s]     = both_called[s] & (chunk[s,i] == chunk[s,j])
    #
    # We use broadcasting: expand dims along axis=2 for i and axis=1 for j
    # → shapes (S, N, 1) and (S, 1, N) → (S, N, N)
    # Sum over the SNP axis (axis=0) → (N, N) matrices.
    # Memory for the intermediate (S, N, N) tensor would be huge for large N,
    # so instead we iterate over i explicitly (N iterations, each O(S*N)).

    N = chunk.shape[1]
    for i in range(N):
        # valid_i : (S,)   calls_i : (S,)
        valid_i = valid[:, i]
        calls_i = chunk[:, i]

        # both_called : (S, N)  — SNPs where sample i AND each sample j are called
        both_called = valid_i[:, np.newaxis] & valid    # (S,1) & (S,N) → (S,N)

        # matched : (S, N)
        matched = both_called & (calls_i[:, np.newaxis] == chunk)

        # Sum over SNPs → (N,) vectors; accumulate into row i
        n_called[i] += both_called.sum(axis=0)
        n_match[i]  += matched.sum(axis=0)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    traw_path = Path(args.traw)
    out_path  = Path(args.out)
    log_path  = out_path.with_suffix(".log")

    if not traw_path.exists():
        sys.exit(f"ERROR: {traw_path} not found")

    print(f"[pairwise_concordance] Input  : {traw_path}", flush=True)
    print(f"[pairwise_concordance] Output : {out_path}",  flush=True)
    print(f"[pairwise_concordance] Chunk  : {args.chunk:,} SNPs", flush=True)

    # ── Parse header ─────────────────────────────────────────────────────────
    sample_ids = read_traw_header(traw_path)
    N = len(sample_ids)
    print(f"[pairwise_concordance] Samples: {N}", flush=True)
    print(f"[pairwise_concordance] Pairs  : {N*(N-1)//2:,}", flush=True)

    # ── Allocate accumulators (int64 to handle ~2M SNPs × 1000 samples) ──────
    n_called = np.zeros((N, N), dtype=np.int64)
    n_match  = np.zeros((N, N), dtype=np.int64)

    # ── Stream and accumulate ─────────────────────────────────────────────────
    t0          = time.time()
    snps_done   = 0
    chunks_done = 0

    for chunk in iter_chunks(traw_path, args.chunk, N):
        accumulate(chunk, n_called, n_match)
        snps_done   += chunk.shape[0]
        chunks_done += 1
        elapsed = time.time() - t0
        print(f"[pairwise_concordance]   chunk {chunks_done:4d} | "
              f"{snps_done:>10,} SNPs | {elapsed:6.1f}s elapsed", flush=True)

    print(f"[pairwise_concordance] Finished streaming {snps_done:,} SNPs "
          f"in {time.time()-t0:.1f}s", flush=True)

    # ── Compute concordance for upper triangle ────────────────────────────────
    rows_i, rows_j, called_ij, match_ij, concordance = [], [], [], [], []

    for i in range(N):
        for j in range(i + 1, N):
            c = int(n_called[i, j])
            m = int(n_match[i, j])
            pct = (m / c * 100.0) if c > 0 else float("nan")
            rows_i.append(sample_ids[i])
            rows_j.append(sample_ids[j])
            called_ij.append(c)
            match_ij.append(m)
            concordance.append(round(pct, 4) if not np.isnan(pct) else None)

    results = pd.DataFrame({
        "SAMPLE_A":        rows_i,
        "SAMPLE_B":        rows_j,
        "N_CALLED":        called_ij,
        "N_MATCH":         match_ij,
        "CONCORDANCE_PCT": concordance,
    })

    results.to_csv(out_path, sep="\t", index=False, na_rep="NA")
    print(f"[pairwise_concordance] Results written → {out_path}", flush=True)

    # ── Log summary ───────────────────────────────────────────────────────────
    valid_pairs = results["CONCORDANCE_PCT"].notna()
    conc        = results.loc[valid_pairs, "CONCORDANCE_PCT"]
    flagged     = results[results["CONCORDANCE_PCT"] >= args.dup_threshold]

    with open(log_path, "w") as log:
        log.write("=" * 60 + "\n")
        log.write("Pairwise Genotype Concordance Summary\n")
        log.write("=" * 60 + "\n")
        log.write(f"Input            : {traw_path}\n")
        log.write(f"Samples          : {N}\n")
        log.write(f"Total pairs      : {len(results):,}\n")
        log.write(f"SNPs processed   : {snps_done:,}\n")
        log.write(f"Chunk size       : {args.chunk:,}\n")
        log.write(f"Dup threshold    : {args.dup_threshold}%\n")
        log.write(f"\nConcordance statistics (all valid pairs):\n")
        log.write(f"  Mean   : {conc.mean():.2f}%\n")
        log.write(f"  Median : {conc.median():.2f}%\n")
        log.write(f"  Min    : {conc.min():.2f}%\n")
        log.write(f"  Max    : {conc.max():.2f}%\n")
        log.write(f"\nPairs >= {args.dup_threshold}% (likely duplicates): {len(flagged)}\n")
        if len(flagged) > 0:
            log.write("\n  SAMPLE_A\tSAMPLE_B\tN_CALLED\tN_MATCH\tCONCORDANCE_PCT\n")
            for _, r in flagged.sort_values("CONCORDANCE_PCT", ascending=False).iterrows():
                log.write(f"  {r.SAMPLE_A}\t{r.SAMPLE_B}\t{r.N_CALLED}\t"
                          f"{r.N_MATCH}\t{r.CONCORDANCE_PCT}\n")
        log.write("\n" + "=" * 60 + "\n")

    print(f"[pairwise_concordance] Log written      → {log_path}", flush=True)
    print(f"[pairwise_concordance] Pairs >= {args.dup_threshold}%: {len(flagged)}", flush=True)

    if len(flagged) > 0:
        print(f"[pairwise_concordance] *** {len(flagged)} likely duplicate pair(s) detected ***",
              flush=True)


if __name__ == "__main__":
    main()