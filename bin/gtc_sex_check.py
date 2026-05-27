#!/usr/bin/env python3
"""
gtc_sex_check.py
================
Compare collected sex (samplesheet / sex_info) with Illumina GTC computed_gender
from gtc_qc_summary.tsv (bcftools +gtc2vcf --extra stats).

Writes outputs compatible with the HTML report and prior CHECK_SEX publish layout.
"""
from __future__ import annotations

import argparse
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

GENDER_MAP = {
    "M": "Male",
    "MALE": "Male",
    "1": "Male",
    "F": "Female",
    "FEMALE": "Female",
    "0": "Female",
    "2": "Female",
}

PLINK_SNPSEX = {"Male": 1, "Female": 2, "Unknown": 0}
SEX_PALETTE = {"Male": "#2196F3", "Female": "#E91E63", "Unknown": "#9E9E9E"}
PLATE_CANDIDATES = [
    "plate", "plate_number", "plate number", "batch", "sentrix_barcode",
    "beadchip_barcode", "barcode",
]


def parse_args():
    p = argparse.ArgumentParser(description="GTC computed_gender vs collected sex")
    p.add_argument("--gtc_qc_summary", required=True)
    p.add_argument("--sex_info", required=True)
    p.add_argument("--samplesheet", default=None,
                   help="Resolved samplesheet CSV (plate / well metadata)")
    p.add_argument("--out_annot", required=True)
    p.add_argument("--out_discord", required=True)
    p.add_argument("--out_multimind", required=True)
    p.add_argument("--out_plate_report", required=True)
    p.add_argument("--out_xy", required=True)
    p.add_argument("--out_plot", required=True)
    p.add_argument("--alert_pct", type=float, default=30.0)
    return p.parse_args()


def normalise_gender(val) -> str:
    if pd.isna(val):
        return "Unknown"
    key = str(val).strip().upper()
    return GENDER_MAP.get(key, "Unknown")


def load_gtc_summary(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = df.columns.str.strip().str.lower()
    if "sample_id" in df.columns:
        df["IID"] = df["sample_id"].astype(str).str.strip()
    elif "iid" in df.columns:
        df["IID"] = df["iid"].astype(str).str.strip()
    else:
        sys.exit("[gtc_sex_check] ERROR: gtc_qc_summary missing sample_id column")
    if "computed_gender" not in df.columns:
        sys.exit("[gtc_sex_check] ERROR: gtc_qc_summary missing computed_gender column")
    df["INFERRED_SEX"] = df["computed_gender"].map(normalise_gender)
    for col in ("p50_x", "p50_y", "logr_deviation", "call_rate"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def load_sex_info(path: str) -> pd.DataFrame:
    sex_info = pd.read_csv(path, sep="\t", dtype=str)
    sex_info.columns = sex_info.columns.str.strip().str.lower()
    sex_info = sex_info.rename(columns={"sampleid": "iid"})
    sex_info["IID"] = sex_info["iid"].astype(str).str.strip()
    sex_info["COLLECTED_SEX"] = sex_info["sex"].map(
        {"0": "Female", "1": "Male"}
    ).fillna("Unknown")
    return sex_info


def load_plate_map(samplesheet_path: str | None) -> pd.DataFrame | None:
    if not samplesheet_path:
        return None
    try:
        ss = pd.read_csv(samplesheet_path, dtype=str)
    except Exception as exc:
        print(f"[gtc_sex_check] WARNING: could not read samplesheet: {exc}", flush=True)
        return None
    ss.columns = ss.columns.str.strip().str.lower()
    id_col = next((c for c in ("sample_id", "sampleid", "iid") if c in ss.columns), None)
    if not id_col:
        return None
    ss["IID"] = ss[id_col].astype(str).str.strip()
    plate_col = next((c for c in PLATE_CANDIDATES if c in ss.columns), None)
    if not plate_col:
        return None
    return ss[["IID", plate_col]].rename(columns={plate_col: "PLATE"})


def is_discordant(row) -> bool:
    if row["INFERRED_SEX"] == "Unknown" or row["COLLECTED_SEX"] == "Unknown":
        return False
    return row["INFERRED_SEX"] != row["COLLECTED_SEX"]


def overall_class(row) -> str:
    if row["INFERRED_SEX"] == "Unknown" or row["COLLECTED_SEX"] == "Unknown":
        return "UNKNOWN"
    if row["DISCORDANT"]:
        return "HARD_DISCORDANT"
    return "CONCORDANT"


def write_plate_report(merged: pd.DataFrame, plate_map: pd.DataFrame | None, out_path: str, alert_pct: float):
    sexcheck = merged.copy()
    if plate_map is not None:
        sexcheck = sexcheck.merge(plate_map, on="IID", how="left")
    else:
        sexcheck["PLATE"] = "UNKNOWN"
    sexcheck["PLATE"] = sexcheck["PLATE"].fillna("UNKNOWN")

    def summarise(grp):
        n_total = len(grp)
        discordant = grp["DISCORDANT"].sum()
        pct = 100 * discordant / n_total if n_total else 0
        iids = grp.loc[grp["DISCORDANT"], "IID"].tolist()
        return pd.Series({
            "N_SAMPLES": n_total,
            "N_DISCORDANT": int(discordant),
            "PCT_DISCORDANT": round(pct, 1),
            "IIDs_DISCORDANT": ",".join(iids),
        })

    report = sexcheck.groupby("PLATE").apply(summarise).reset_index()
    report = report.sort_values("PCT_DISCORDANT", ascending=False)
    report.to_csv(out_path, sep="\t", index=False)

    print("[gtc_sex_check] Per-plate discordance:", flush=True)
    for _, row in report.iterrows():
        flag = " *** POSSIBLE PLATE SWAP ***" if row["PCT_DISCORDANT"] >= alert_pct else ""
        print(
            f"  Plate {row['PLATE']}: {int(row['N_DISCORDANT'])}/{int(row['N_SAMPLES'])} "
            f"discordant ({row['PCT_DISCORDANT']:.1f}%){flag}",
            flush=True,
        )


def write_xy_tsv(merged: pd.DataFrame, out_path: str):
    cols = ["IID", "COLLECTED_SEX", "INFERRED_SEX"]
    xy = merged[cols].copy()
    if "p50_x" in merged.columns:
        xy["MEAN_X_X"] = merged["p50_x"]
    if "p50_y" in merged.columns:
        xy["MEAN_Y_Y"] = merged["p50_y"]
    xy.to_csv(out_path, sep="\t", index=False)


def write_plot(merged: pd.DataFrame, out_path: str):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    fig.suptitle("Sex Check – GTC computed gender", fontsize=12, fontweight="bold")

    ax = axes[0]
    if "logr_deviation" in merged.columns:
        for sex, grp in merged.groupby("COLLECTED_SEX"):
            ax.hist(
                grp["logr_deviation"].dropna(),
                bins=40,
                alpha=0.65,
                color=SEX_PALETTE.get(sex, "#9E9E9E"),
                label=sex,
                edgecolor="white",
                linewidth=0.3,
            )
        ax.set_xlabel("log R deviation (GTC)", fontsize=10)
        ax.set_ylabel("Samples", fontsize=10)
        ax.set_title("logR deviation by collected sex", fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(True, linewidth=0.4, alpha=0.5)
    else:
        ax.set_visible(False)

    ax = axes[1]
    if "p50_x" in merged.columns and "p50_y" in merged.columns:
        for sex, grp in merged.groupby("COLLECTED_SEX"):
            ax.scatter(
                grp["p50_x"],
                grp["p50_y"],
                c=SEX_PALETTE.get(sex, "#9E9E9E"),
                alpha=0.55,
                s=12,
                edgecolors="none",
                label=sex,
            )
        ax.set_xlabel("Median X intensity (p50_x)", fontsize=10)
        ax.set_ylabel("Median Y intensity (p50_y)", fontsize=10)
        ax.set_title("X/Y intensity by collected sex", fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(True, linewidth=0.4, alpha=0.5)
    else:
        ax.set_visible(False)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()

    gtc = load_gtc_summary(args.gtc_qc_summary)
    sex_info = load_sex_info(args.sex_info)
    plate_map = load_plate_map(args.samplesheet)

    merged = gtc.merge(sex_info[["IID", "COLLECTED_SEX"]], on="IID", how="left")
    merged["COLLECTED_SEX"] = merged["COLLECTED_SEX"].fillna("Unknown")
    merged["DISCORDANT"] = merged.apply(is_discordant, axis=1)
    merged["OVERALL_CLASS"] = merged.apply(overall_class, axis=1)

    merged["FID"] = merged["IID"]
    merged["PEDSEX"] = merged["COLLECTED_SEX"].map({"Male": 1, "Female": 2}).fillna(0).astype(int)
    merged["SNPSEX"] = merged["INFERRED_SEX"].map(PLINK_SNPSEX).fillna(0).astype(int)
    merged["STATUS"] = merged["DISCORDANT"].map({True: "PROBLEM", False: "OK"})
    merged.loc[merged["INFERRED_SEX"] == "Unknown", "STATUS"] = "PROBLEM"
    merged["F"] = merged["logr_deviation"] if "logr_deviation" in merged.columns else float("nan")
    merged["INFERENCE_METHOD"] = "GTC_computed_gender"

    annot_cols = [
        "FID", "IID", "PEDSEX", "SNPSEX", "STATUS", "F",
        "COLLECTED_SEX", "INFERRED_SEX", "DISCORDANT", "INFERENCE_METHOD",
        "computed_gender",
    ]
    annot_cols = [c for c in annot_cols if c in merged.columns]
    merged[annot_cols].to_csv(args.out_annot, sep="\t", index=False)

    discord = merged[merged["DISCORDANT"]][
        ["FID", "IID", "COLLECTED_SEX", "INFERRED_SEX", "F", "STATUS"]
    ]
    discord.to_csv(args.out_discord, sep="\t", index=False)

    multimind = merged[["IID", "OVERALL_CLASS", "F"]].rename(columns={"F": "F_base"})
    multimind.to_csv(args.out_multimind, sep="\t", index=False)

    write_plate_report(merged, plate_map, args.out_plate_report, args.alert_pct)
    write_xy_tsv(merged, args.out_xy)
    write_plot(merged, args.out_plot)

    n_missing_gtc = len(sex_info) - merged["computed_gender"].notna().sum() if "computed_gender" in merged.columns else 0
    if len(sex_info) > len(merged):
        n_missing_gtc = len(sex_info) - len(merged)

    print(f"[gtc_sex_check] Samples in GTC summary : {len(gtc)}", flush=True)
    print(f"[gtc_sex_check] Samples with sex_info  : {len(sex_info)}", flush=True)
    print(f"[gtc_sex_check] Sex-discordant        : {int(merged['DISCORDANT'].sum())}", flush=True)
    print(f"[gtc_sex_check] STATUS=PROBLEM          : {(merged['STATUS'] == 'PROBLEM').sum()}", flush=True)
    if n_missing_gtc:
        print(f"[gtc_sex_check] WARNING: {n_missing_gtc} sex_info sample(s) missing from GTC summary", flush=True)


if __name__ == "__main__":
    main()
