#!/usr/bin/env python3
"""
filter_chrx_snps.py
H3AGWAS stats_x.py-style sex-stratified chrX SNP QC.
Reads PLINK --freq / --missing / --hardy outputs and writes SNPs passing all filters.
"""
import argparse
import sys

import pandas as pd
import scipy.stats as stats

PLINK_SEP = r"\s+"


def parse_args():
    parser = argparse.ArgumentParser(description="Filter chrX SNPs for sex check")
    parser.add_argument("--male_prefix", required=True,
                        help="PLINK output prefix for male freq/miss (e.g. chrX_male)")
    parser.add_argument("--female_prefix", required=True,
                        help="PLINK output prefix for female freq/miss/HWE (e.g. chrX_female)")
    parser.add_argument("--bim", required=True,
                        help="chrX_hhmissing.bim used for fallback SNP list")
    parser.add_argument("--out_summary", required=True)
    parser.add_argument("--out_snps", required=True)
    parser.add_argument("--x_maf_male", type=float, default=0.01)
    parser.add_argument("--x_maf_female", type=float, default=0.01)
    parser.add_argument("--x_miss_male", type=float, default=0.8)
    parser.add_argument("--x_miss_female", type=float, default=0.05)
    parser.add_argument("--x_diff_miss", type=float, default=0.6)
    parser.add_argument("--x_hwe_female", type=float, default=1e-4)
    parser.add_argument("--x_fisher_p", type=float, default=1e-4)
    return parser.parse_args()


def load_plink_table(prefix, suffix):
    return pd.read_csv(f"{prefix}.{suffix}", sep=PLINK_SEP, engine="python")


def write_fallback_snps(bim_path, out_snps):
    bim = pd.read_csv(
        bim_path, sep="\t", header=None,
        names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
    )
    bim[["SNP"]].to_csv(out_snps, index=False, header=False)
    print(f"[chrX_QC] Fallback: wrote all {len(bim)} variants to {out_snps}", flush=True)


def build_joint_dataset(male_prefix, female_prefix):
    frq_m = load_plink_table(male_prefix, "frq")
    lmiss_m = load_plink_table(male_prefix, "lmiss")
    frq_f = load_plink_table(female_prefix, "frq")
    lmiss_f = load_plink_table(female_prefix, "lmiss")

    frq_m = frq_m.rename(columns={
        "MAF": "MAF_male", "NCHROBS": "NCHROBS_male", "A1": "A1_male", "A2": "A2_male",
    })
    frq_m["N_GENO_male"] = frq_m.get("N_GENO", frq_m["NCHROBS_male"])
    lmiss_m = lmiss_m.rename(columns={"F_MISS": "MISS_male"})

    data = frq_m[["CHR", "SNP", "MAF_male", "NCHROBS_male", "N_GENO_male", "A1_male", "A2_male"]].merge(
        lmiss_m[["CHR", "SNP", "MISS_male"]], on=["CHR", "SNP"],
    )

    frq_f = frq_f.rename(columns={
        "MAF": "MAF_female", "NCHROBS": "NCHROBS_female", "A1": "A1_female", "A2": "A2_female",
    })
    frq_f["N_GENO_female"] = frq_f.get(
        "N_GENO", (frq_f["NCHROBS_female"] / 2).round().astype(int),
    )
    lmiss_f = lmiss_f.rename(columns={"F_MISS": "MISS_female"})

    data = data.merge(
        frq_f[["CHR", "SNP", "MAF_female", "N_GENO_female", "NCHROBS_female", "A1_female", "A2_female"]],
        on=["CHR", "SNP"],
    ).merge(
        lmiss_f[["CHR", "SNP", "MISS_female"]], on=["CHR", "SNP"],
    )
    return data


def attach_hwe(data, female_prefix):
    try:
        hwe = load_plink_table(female_prefix, "hwe")
        if "TEST" in hwe.columns:
            hwe = hwe[hwe["TEST"] == "ALL"]
        hwe_cols = ["CHR", "SNP"] + (["P"] if "P" in hwe.columns else [])
        data = data.merge(hwe[hwe_cols], on=["CHR", "SNP"], how="left")
        data = data.rename(columns={"P": "HWE_P_female"})
        print(f"[chrX_QC] Female HWE loaded: {len(hwe)} SNPs", flush=True)
    except Exception as exc:
        print(f"[chrX_QC] WARNING: female HWE unavailable: {exc}", flush=True)
        data["HWE_P_female"] = 1.0
    return data


def fisher_p(row):
    vals = [row["N_A1_male"], row["N_A2_male"], row["N_A1_female"], row["N_A2_female"]]
    if any(pd.isna(v) or int(v) < 0 for v in vals):
        return float("nan")
    try:
        _, p = stats.fisher_exact([
            [int(vals[0]), int(vals[1])],
            [int(vals[2]), int(vals[3])],
        ])
        return p
    except Exception:
        return float("nan")


def apply_filters(data, args):
    maf_m, maf_f = args.x_maf_male, args.x_maf_female
    data["pass_maf_male"] = data["MAF_male"].between(maf_m, 1 - maf_m, inclusive="both")
    data["pass_maf_female"] = data["MAF_female"].between(maf_f, 1 - maf_f, inclusive="both")
    data["pass_miss_male"] = data["MISS_male"] < args.x_miss_male
    data["pass_miss_female"] = data["MISS_female"] < args.x_miss_female
    data["pass_hwe_female"] = data["HWE_P_female"].fillna(1.0) >= args.x_hwe_female
    data["pass_diff_miss"] = (data["MISS_male"] - data["MISS_female"]).abs() < args.x_diff_miss

    flip = data["A1_male"] != data["A1_female"]
    data["MAF_male_adj"] = data["MAF_male"].copy()
    data.loc[flip, "MAF_male_adj"] = 1 - data.loc[flip, "MAF_male"]
    data["N_A1_male"] = (data["MAF_male_adj"] * data["N_GENO_male"]).round().astype("Int64")
    data["N_A2_male"] = data["N_GENO_male"] - data["N_A1_male"]
    data["N_A1_female"] = (data["MAF_female"] * data["N_GENO_female"]).round().astype("Int64")
    data["N_A2_female"] = data["N_GENO_female"] - data["N_A1_female"]

    data["fisher_p"] = data.apply(fisher_p, axis=1)
    data["pass_fisher"] = data["fisher_p"].fillna(1.0) > args.x_fisher_p

    filter_cols = [
        "pass_maf_male", "pass_maf_female", "pass_miss_male",
        "pass_miss_female", "pass_hwe_female", "pass_diff_miss", "pass_fisher",
    ]
    data["PASS_ALL"] = data[filter_cols].all(axis=1)
    return data, filter_cols


def main():
    args = parse_args()
    print(
        f"[chrX_QC] Thresholds — MAF_M>={args.x_maf_male}, MAF_F>={args.x_maf_female}, "
        f"MISS_M(post)<{args.x_miss_male}, MISS_F<{args.x_miss_female}, "
        f"|dMISS|<{args.x_diff_miss}, HWE_F>={args.x_hwe_female}, Fisher>{args.x_fisher_p}",
        flush=True,
    )

    try:
        load_plink_table(args.male_prefix, "frq")
        load_plink_table(args.male_prefix, "lmiss")
        load_plink_table(args.female_prefix, "frq")
        load_plink_table(args.female_prefix, "lmiss")
    except Exception as exc:
        print(f"[chrX_QC] WARNING: freq/miss unavailable ({exc}) — fallback to all variants.", flush=True)
        write_fallback_snps(args.bim, args.out_snps)
        sys.exit(0)

    data = build_joint_dataset(args.male_prefix, args.female_prefix)
    print(f"[chrX_QC] Joint M+F dataset: {len(data)} SNPs", flush=True)

    data = attach_hwe(data, args.female_prefix)
    data, filter_cols = apply_filters(data, args)

    n_pass = int(data["PASS_ALL"].sum())
    print("[chrX_QC] Filter summary:", flush=True)
    for col in filter_cols:
        print(f"  Failed {col:<22}: {(~data[col]).sum()}", flush=True)
    print(f"  SNPs passing all filters     : {n_pass} / {len(data)}", flush=True)

    if n_pass == 0:
        print("[chrX_QC] WARNING: 0 SNPs passed — fallback to all variants.", flush=True)
        write_fallback_snps(args.bim, args.out_snps)
        sys.exit(0)

    data.to_csv(args.out_summary, index=False)
    data.loc[data["PASS_ALL"], ["SNP"]].to_csv(args.out_snps, index=False, header=False)
    print(f"[chrX_QC] Written: {args.out_snps}", flush=True)


if __name__ == "__main__":
    main()
