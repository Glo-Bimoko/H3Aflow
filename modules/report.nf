process REPORT {
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path sexcheck
    path multimind
    path plate_report
    path xy_tsv
    path qc_stats
    path genome
    path eigenvec
    path sex_info
    path concordance   // pairwise_concordance.tsv  (pass file("NO_FILE") if absent)
    path samplesheet   // original samplesheet CSV  (pass file("NO_FILE") if absent)
    // GTC-level QC summary and poor-GC list (optional)
    path gtc_qc_summary
    path poorgc10

    output:
    path "qc_report.html"
    path "flagged_samples.tsv"
    path "cohort_summary.tsv"
    path "sexcheck_multimind.tsv"
    path "sexcheck_plate_report.tsv"

    script:
    def conc_arg  = concordance.name != "NO_FILE" ? "--concordance ${concordance}" : ""
    def sheet_arg = samplesheet.name != "NO_FILE"  ? "--samplesheet ${samplesheet}" : ""
    def gtc_arg   = gtc_qc_summary.name != "NO_FILE" ? "--gtc_qc_summary ${gtc_qc_summary}" : ""
    def poor_arg  = poorgc10.name != "NO_FILE" ? "--poorgc10 ${poorgc10}" : ""
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
        ${sheet_arg} \\
        ${gtc_arg}  \\
        ${poor_arg} \\
        --gc10_threshold ${params.gc10_threshold ?: 0.15}
    """
}