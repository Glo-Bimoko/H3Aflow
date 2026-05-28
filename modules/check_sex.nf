process CHECK_SEX {
    tag "cohort"
    publishDir "${params.outdir}/qc/check_sex", mode: 'copy'

    input:
    path gtc_qc_summary
    path sex_info
    path resolved_samplesheet

    output:
    path "sexcheck.txt",              emit: sexcheck
    path "sex_discordance.txt",       emit: discordance
    path "sex_plot.png",              emit: plot
    path "sex_check.log",             emit: log
    path "sexcheck_multimind.tsv",    emit: multimind
    path "sexcheck_plate_report.tsv", emit: plate_report
    path "xy_intensity.tsv",          emit: xy_tsv

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG=sex_check.log
    exec > >(tee -a "\$LOG") 2>&1

    echo "========================================"
    echo "Sex Check (GTC computed_gender): \$(date)"
    echo "========================================"

    ${params.python} ${projectDir}/bin/gtc_sex_check.py \\
        --gtc_qc_summary   ${gtc_qc_summary} \\
        --sex_info         ${sex_info} \\
        --samplesheet      ${resolved_samplesheet} \\
        --out_annot        sexcheck.txt \\
        --out_discord      sex_discordance.txt \\
        --out_multimind    sexcheck_multimind.tsv \\
        --out_plate_report sexcheck_plate_report.tsv \\
        --out_xy           xy_intensity.tsv \\
        --out_plot         sex_plot.png

    echo "========================================"
    echo "Sex Check Completed: \$(date)"
    echo "========================================"
    """
}