process REPORT {
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path sexcheck
    path xy_tsv
    path qc_stats
    path genome
    path eigenvec
    path sex_info
    path concordance   // pairwise_concordance.tsv  (pass file("NO_FILE") if absent)
    path samplesheet   // original samplesheet CSV  (pass file("NO_FILE") if absent)

    output:
    path "qc_report.html"
    path "flagged_samples.tsv"
    path "cohort_summary.tsv"

    script:
    def conc_arg  = concordance.name != "NO_FILE" ? "--concordance ${concordance}" : ""
    def sheet_arg = samplesheet.name != "NO_FILE"  ? "--samplesheet ${samplesheet}" : ""
    """
    python ${projectDir}/bin/generate_report.py \\
        --sexcheck    ${sexcheck}   \\
        --xy_tsv      ${xy_tsv}     \\
        --qc_stats    ${qc_stats}   \\
        --genome      ${genome}     \\
        --eigenvec    ${eigenvec}   \\
        --sex_info    ${sex_info}   \\
        --out_html    qc_report.html \\
        --out_flagged flagged_samples.tsv \\
        --out_summary cohort_summary.tsv \\
        --concordance_warn  ${params.concordance_warn  ?: 84}  \\
        --concordance_flag  ${params.concordance_flag  ?: 99}  \\
        ${conc_arg}  \\
        ${sheet_arg}
    """
}