process REPORT {
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path sexcheck
    path xy_tsv
    path qc_stats
    path genome
    path eigenvec
    path sex_info

    output:
    path "qc_report.html"
    path "flagged_samples.tsv"
    path "cohort_summary.tsv"

    script:
    """
    python ${projectDir}/bin/generate_report.py \\
        --sexcheck  ${sexcheck} \\
        --xy_tsv    ${xy_tsv} \\
        --qc_stats  ${qc_stats} \\
        --genome    ${genome} \\
        --eigenvec  ${eigenvec} \\
        --sex_info  ${sex_info} \\
        --out_html  qc_report.html \\
        --out_flagged flagged_samples.tsv \\
        --out_summary cohort_summary.tsv
    """
}