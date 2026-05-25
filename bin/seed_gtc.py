#!/usr/bin/env python3
"""
seed_gtc.py
Copy a pre-seeded per-sample GTC into the Nextflow work directory.
Used for intentional GTC-level duplicates (same .gtc, different Sample ID).
"""
import argparse
import shutil
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Stage a pre-seeded GTC for one sample")
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--gtc_dir", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    source = Path(args.gtc_dir) / f"{args.sample_id}.gtc"
    output = Path(args.output)

    if not source.is_file():
        sys.exit(f"[seed_gtc] ERROR: pre-seeded GTC not found: {source}")

    if source.resolve() != output.resolve():
        shutil.copy2(source, output)

    print(f"[seed_gtc] Staged seeded GTC for {args.sample_id} from {source}", flush=True)


if __name__ == "__main__":
    main()
