#!/usr/bin/env python3
"""
sexcheck_multimind.py
H3AGWAS xCheck.py-style multi-missingness sex check.
Re-runs PLINK --check-sex at multiple --mind thresholds and classifies each sample.
"""
import argparse
import os
import subprocess
import sys

import pandas as pd

PLINK_SEP = r"\s+"
ID_COLS = {"FID": object, "IID": object}


def parse_args():
    parser = argparse.ArgumentParser(description="Multi-threshold PLINK sex check")
    parser.add_argument("--bfile", required=True, help="PLINK bfile prefix (no extension)")
    parser.add_argument("--out", required=True, help="Output TSV path")
    parser.add_argument("--f_hi_female", type=float, default=0.2)
    parser.add_argument("--f_lo_male", type=float, default=0.8)
    parser.add_argument("--mind_thresholds", default="0.01,0.03,0.05,0.1",
                        help="Comma-separated --mind values")
    return parser.parse_args()


def parse_thresholds(raw):
    return [float(x.strip()) for x in raw.split(",") if x.strip()]


def run_checksex(bfile, f_hi_female, f_lo_male, mind, out_prefix):
    cmd = [
        "plink", "--bfile", bfile,
        "--check-sex", str(f_hi_female), str(f_lo_male),
        "--out", out_prefix,
    ]
    if mind is not None:
        cmd += ["--mind", str(mind)]
    subprocess.run(cmd, capture_output=True, text=True, check=False)

    outfile = f"{out_prefix}.sexcheck"
    if not os.path.exists(outfile):
        print(f"  [mind={mind}] PLINK produced no output — all samples excluded?", flush=True)
        return None
    return pd.read_csv(outfile, sep=PLINK_SEP, dtype=ID_COLS)


def classify(row, f_hi_female, f_lo_male):
    if row["STATUS"] == "OK":
        return "OK"
    f_val = row["F"]
    if pd.isna(f_val):
        return "S"
    if f_hi_female < f_val < f_lo_male:
        return "S"
    return "H"


def overall_class(row, mind_cols):
    vals = [row[c] for c in mind_cols if row[c] != "X"]
    if not vals:
        return "EXCLUDED"
    if all(v == "OK" for v in vals):
        return "CONCORDANT"
    if all(v == "H" for v in vals):
        return "HARD_DISCORDANT"
    if any(v == "H" for v in vals):
        return "DISCORDANT_VARIABLE"
    return "SUSPECT"


def print_classification_summary(result, mind_cols):
    print("", flush=True)
    print("[multi_mind] Overall classification summary:", flush=True)
    for cls, cnt in result["OVERALL_CLASS"].value_counts().items():
        print(f"  {cls:<25}: {cnt}", flush=True)

    print("", flush=True)
    print("[multi_mind] Hard discordant samples (wrong sex at ALL thresholds):", flush=True)
    hard = result[result["OVERALL_CLASS"] == "HARD_DISCORDANT"]
    if hard.empty:
        print("  None — all discordant samples are missingness-driven.", flush=True)
    else:
        for _, row in hard.iterrows():
            print(
                f"  IID={row['IID']}  PEDSEX={row['PEDSEX']}  F_base={row['F_base']:.4f}",
                flush=True,
            )

    print("", flush=True)
    print("[multi_mind] Variable discordant (missingness-driven):", flush=True)
    variable = result[result["OVERALL_CLASS"] == "DISCORDANT_VARIABLE"]
    if variable.empty:
        print("  None.", flush=True)
    else:
        for _, row in variable.iterrows():
            vals_str = "  ".join(f"{col}={row[col]}" for col in mind_cols)
            print(f"  IID={row['IID']}  F_base={row['F_base']:.4f}  {vals_str}", flush=True)


def main():
    args = parse_args()
    thresholds = parse_thresholds(args.mind_thresholds)

    print(f"[multi_mind] Running --check-sex at --mind thresholds: {thresholds}", flush=True)
    print(
        f"[multi_mind] F thresholds: male>{args.f_lo_male}, female<{args.f_hi_female}",
        flush=True,
    )

    base_df = run_checksex(args.bfile, args.f_hi_female, args.f_lo_male, None, "sexcheck_mm_base")
    if base_df is None:
        print("[multi_mind] ERROR: base run failed.", flush=True)
        pd.DataFrame(columns=["FID", "IID", "mind_base"]).to_csv(args.out, sep="\t", index=False)
        sys.exit(0)

    base_df["class_base"] = base_df.apply(
        classify, axis=1, args=(args.f_hi_female, args.f_lo_male),
    )
    result = base_df[["FID", "IID", "PEDSEX", "F", "STATUS", "class_base"]].copy()
    result = result.rename(columns={"class_base": "mind_base", "F": "F_base"})

    for mind in thresholds:
        label = f"mind_{mind}"
        df = run_checksex(args.bfile, args.f_hi_female, args.f_lo_male, mind, f"sexcheck_mm_{label}")
        if df is None:
            result[label] = "X"
            print(f"  [mind={mind}] All samples excluded by --mind.", flush=True)
            continue

        df["cls"] = df.apply(classify, axis=1, args=(args.f_hi_female, args.f_lo_male))
        cls_map = df.set_index("IID")["cls"].to_dict()
        result[label] = result["IID"].map(cls_map).fillna("X")

        n_ok = (result[label] == "OK").sum()
        n_s = (result[label] == "S").sum()
        n_h = (result[label] == "H").sum()
        n_x = (result[label] == "X").sum()
        print(
            f"  [mind={mind}] {len(df)} samples kept — OK:{n_ok}  S:{n_s}  H:{n_h}  X:{n_x}",
            flush=True,
        )

    mind_cols = ["mind_base"] + [f"mind_{m}" for m in thresholds]
    result["OVERALL_CLASS"] = result.apply(overall_class, axis=1, args=(mind_cols,))

    print_classification_summary(result, mind_cols)
    result.to_csv(args.out, sep="\t", index=False)
    print("", flush=True)
    print(f"[multi_mind] Written: {args.out}", flush=True)


if __name__ == "__main__":
    main()
