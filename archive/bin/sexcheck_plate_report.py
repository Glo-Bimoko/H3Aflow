#!/usr/bin/env python3
"""
sexcheck_plate_report.py
Summarises sex-check discordance rates per plate to detect batch-level label swaps.
"""
import argparse
import sys

import pandas as pd

PLINK_SEP = r"\s+"
DISCORDANT_CLASSES = ["HARD_DISCORDANT", "DISCORDANT_VARIABLE", "SUSPECT"]
PLATE_CANDIDATES = [
    "plate", "plate_number", "batch", "sentrix_barcode",
    "beadchip_barcode", "barcode",
]


def parse_args():
    parser = argparse.ArgumentParser(description="Per-plate sex-check discordance report")
    parser.add_argument("--sexcheck", required=True, help="PLINK sexcheck.sexcheck file")
    parser.add_argument("--multimind", required=True, help="sexcheck_multimind.tsv")
    parser.add_argument("--sex_info", required=True, help="sex_info.tsv with sample metadata")
    parser.add_argument("--out", required=True)
    parser.add_argument("--alert_pct", type=float, default=30.0,
                        help="Discordance % threshold for plate-swap alert")
    return parser.parse_args()


def load_sexcheck(path):
    return pd.read_csv(path, sep=PLINK_SEP, dtype={"FID": object, "IID": object})


def attach_overall_class(sexcheck, multimind_path):
    try:
        multimind = pd.read_csv(multimind_path, sep="\t", dtype={"FID": object, "IID": object})
        return sexcheck.merge(multimind[["IID", "OVERALL_CLASS", "F_base"]], on="IID", how="left")
    except Exception:
        sexcheck = sexcheck.copy()
        sexcheck["OVERALL_CLASS"] = sexcheck["STATUS"].map({
            "PROBLEM": "HARD_DISCORDANT",
            "OK": "CONCORDANT",
        })
        return sexcheck


def find_plate_column(sex_info):
    sex_info = sex_info.copy()
    sex_info.columns = [c.lower().strip() for c in sex_info.columns]
    sex_info = sex_info.rename(columns={"sampleid": "IID"})
    sex_info["IID"] = sex_info["IID"].astype(str).str.strip()
    for candidate in PLATE_CANDIDATES:
        if candidate in sex_info.columns:
            return sex_info, candidate
    return sex_info, None


def attach_plate(sexcheck, sex_info_path):
    try:
        sex_info, plate_col = find_plate_column(pd.read_csv(sex_info_path, sep="\t", dtype=str))
        if plate_col:
            merged = sexcheck.merge(sex_info[["IID", plate_col]], on="IID", how="left")
            merged = merged.rename(columns={plate_col: "PLATE"})
            merged["PLATE"] = merged["PLATE"].fillna("UNKNOWN")
            return merged

        print(
            "[plate_report] No plate column found in sex_info — using IID numeric prefix as proxy.",
            flush=True,
        )
    except Exception as exc:
        print(f"[plate_report] Could not load sex_info for plate join: {exc}", flush=True)

    sexcheck = sexcheck.copy()
    sexcheck["PLATE"] = sexcheck["IID"].astype(str).str[:4].fillna("UNKNOWN")
    return sexcheck


def summarise_plate(grp):
    n_total = len(grp)
    discordant = grp["OVERALL_CLASS"].isin(DISCORDANT_CLASSES).sum()
    pct = 100 * discordant / n_total if n_total > 0 else 0
    iids = grp.loc[grp["OVERALL_CLASS"].isin(DISCORDANT_CLASSES), "IID"].tolist()
    return pd.Series({
        "N_SAMPLES": n_total,
        "N_DISCORDANT": discordant,
        "PCT_DISCORDANT": round(pct, 1),
        "IIDs_DISCORDANT": ",".join(iids),
    })


def print_report(report, alert_pct):
    print("[plate_report] Per-plate discordance:", flush=True)
    for _, row in report.iterrows():
        flag = " *** POSSIBLE PLATE SWAP ***" if row["PCT_DISCORDANT"] >= alert_pct else ""
        print(
            f"  Plate {row['PLATE']}: {int(row['N_DISCORDANT'])}/{int(row['N_SAMPLES'])} "
            f"discordant ({row['PCT_DISCORDANT']:.1f}%){flag}",
            flush=True,
        )

    high_discord = report[report["PCT_DISCORDANT"] >= alert_pct]
    if not high_discord.empty:
        print("", flush=True)
        print(f"[plate_report] ALERT: {len(high_discord)} plate(s) with >={alert_pct:.0f}% discordance.", flush=True)
        print("[plate_report] This strongly suggests a plate-level sex label swap.", flush=True)
        print("[plate_report] Cross-check these plates against the original samplesheet.", flush=True)
    else:
        print(f"[plate_report] No plates with >={alert_pct:.0f}% discordance (no obvious plate swap).", flush=True)


def write_empty_report(out_path):
    pd.DataFrame(columns=["PLATE", "N_SAMPLES", "N_DISCORDANT", "PCT_DISCORDANT"]).to_csv(
        out_path, sep="\t", index=False,
    )


def main():
    args = parse_args()

    try:
        sexcheck = load_sexcheck(args.sexcheck)
    except Exception as exc:
        print(f"[plate_report] Could not load sexcheck: {exc}", flush=True)
        write_empty_report(args.out)
        sys.exit(0)

    sexcheck = attach_overall_class(sexcheck, args.multimind)
    sexcheck = attach_plate(sexcheck, args.sex_info)

    report = sexcheck.groupby("PLATE").apply(summarise_plate).reset_index()
    report = report.sort_values("PCT_DISCORDANT", ascending=False)
    report.to_csv(args.out, sep="\t", index=False)

    print_report(report, args.alert_pct)


if __name__ == "__main__":
    main()
