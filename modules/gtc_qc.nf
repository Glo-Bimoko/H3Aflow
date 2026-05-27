process GTC_QC {
    tag "${plate}"

    input:
    tuple val(plate), path(gtcs)
    val bpm
    val egt
    val fasta

    output:
    path("${plate}_gtc_stats.tsv"), emit: stats

    script:
    def safe_plate = plate.replaceAll(/[^A-Za-z0-9_\-]/, '_')
    """
    set -euo pipefail
    export BCFTOOLS_PLUGINS="${projectDir}/bin/plugins/gtc2vcf\${BCFTOOLS_PLUGINS:+:\${BCFTOOLS_PLUGINS}}"

    # Same pattern as GTC_TO_VCF: --gtcs takes a list file or directory of many .gtc files.
    ls -1 *.gtc > gtc_list.txt
    echo "[GTC_QC] Plate ${plate}: \$(wc -l < gtc_list.txt) GTC files"

    bcftools +gtc2vcf \\
        --no-version \\
        --bpm       "${bpm}" \\
        --egt       "${egt}" \\
        --fasta-ref "${fasta}" \\
        --gtcs      gtc_list.txt \\
        --extra     "${safe_plate}_gtc_stats.tsv" \\
        -Ou > /dev/null

    if [ "${safe_plate}" != "${plate}" ]; then
        mv "${safe_plate}_gtc_stats.tsv" "${plate}_gtc_stats.tsv"
    fi

    if [ ! -s "${plate}_gtc_stats.tsv" ]; then
        echo "[GTC_QC] ERROR: empty stats for plate ${plate}" >&2
        exit 1
    fi
    """
}

process FILTER_GTC_SAMPLES {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path(stats_files)
    val(gc10_threshold)
    val(call_rate_threshold)

    output:
    path("poorgc10.lst"), emit: poor_gc10_list
    path("gtc_qc_summary.tsv"), emit: summary
    path("gtc_qc_flags.tsv"), emit: flags

    script:
    """
    python3 << 'EOF'
import pandas as pd
import glob
import os
RENAME = {
    "gencall_score_10_percentile": "p10_gc",
    "gencall_score_50_percentile": "p50_gc",
}

def sample_id_from_row(df):
    if "gtc" in df.columns:
        return (
            df["gtc"].astype(str).str.strip()
            .str.replace(".gtc", "", regex=False)
            .str.strip()
        )
    if "sample_id" in df.columns:
        return df["sample_id"].astype(str).str.strip()
    return None

all_stats = []
for f in sorted(glob.glob("*_gtc_stats.tsv")):
    if os.path.getsize(f) == 0:
        continue
    try:
        df = pd.read_csv(f, sep="\\t")
        ids = sample_id_from_row(df)
        if ids is None:
            print(f"Warning: no sample id column in {f}", flush=True)
            continue
        df["sample_id"] = ids
        all_stats.append(df)
    except Exception as e:
        print(f"Warning: Could not parse {f}: {e}", flush=True)

if all_stats:
    stats_df = pd.concat(all_stats, ignore_index=True)
    stats_df = stats_df.rename(columns=RENAME)

    cols_of_interest = [
        "sample_id", "call_rate", "p10_gc", "p50_gc",
        "computed_gender", "logr_deviation",
        "p05_x", "p50_x", "p95_x", "p05_y", "p50_y", "p95_y",
    ]
    cols_available = ["sample_id"] + [c for c in cols_of_interest[1:] if c in stats_df.columns]
    stats_df = stats_df[cols_available]

    if "p10_gc" in stats_df.columns:
        stats_df["pass_gc10"] = stats_df["p10_gc"] >= ${gc10_threshold}
    else:
        stats_df["pass_gc10"] = True

    if "call_rate" in stats_df.columns:
        stats_df["pass_cr"] = stats_df["call_rate"] >= ${call_rate_threshold}
    else:
        stats_df["pass_cr"] = True

    stats_df["pass_qc"] = stats_df["pass_gc10"] & stats_df["pass_cr"]

    stats_df.to_csv("gtc_qc_summary.tsv", sep="\\t", index=False)

    flags_df = stats_df[~stats_df["pass_qc"]][["sample_id", "pass_gc10", "pass_cr"]].copy()
    flags_df.to_csv("gtc_qc_flags.tsv", sep="\\t", index=False)

    failed = stats_df[~stats_df["pass_gc10"]]["sample_id"].values
    with open("poorgc10.lst", "w") as fh:
        for s in failed:
            fh.write(f"{s} {s}\\n")

    print(
        f"GTC QC: {len(stats_df)} samples, {len(failed)} below p10_gc threshold (${gc10_threshold})",
        flush=True,
    )
else:
    open("poorgc10.lst", "w").close()
    open("gtc_qc_summary.tsv", "w").close()
    open("gtc_qc_flags.tsv", "w").close()
    print("Warning: No GTC stats files found or all were empty", flush=True)
EOF
    """
}
