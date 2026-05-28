process GENOTYPE_CONCORDANCE {
    tag "cohort"
    publishDir "${params.outdir}/qc/concordance", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)
    path samplesheet

    output:
    path "pairwise_concordance.tsv",  emit: concordance
    path "pairwise_concordance.log",  emit: log
    path "idat_duplicate_report.tsv", emit: idat_duplicates

    script:
    def prefix        = bed.baseName
    def chunk         = params.concordance_chunk      ?: 50000
    def dup_threshold = params.concordance_dup_thresh ?: 99.5

    """
    echo "========================================" > concordance_run.log
    echo "GENOTYPE_CONCORDANCE Started: \$(date)"  >> concordance_run.log
    echo "========================================"  >> concordance_run.log
    echo "Input prefix : ${prefix}"                 >> concordance_run.log
    echo "Samples      : \$(wc -l < ${fam})"        >> concordance_run.log
    echo "SNPs         : \$(wc -l < ${bim})"        >> concordance_run.log

    plink \\
        --bfile   ${prefix}        \\
        --recode  A-transpose      \\
        --autosome-num 26          \\
        --allow-extra-chr          \\
        --allow-no-sex             \\
        --out     cohort_for_conc  \\
        2>&1 | tee -a concordance_run.log

    ${params.python} ${projectDir}/bin/pairwise_concordance.py \\
        --traw          cohort_for_conc.traw     \\
        --out           pairwise_concordance.tsv \\
        --chunk         ${chunk}                 \\
        --dup-threshold ${dup_threshold}          \\
        2>&1 | tee -a concordance_run.log

    ${params.python} ${projectDir}/bin/check_idat_duplicates.py \\
        --samplesheet ${samplesheet}             \\
        --fam         ${fam}                     \\
        --concordance pairwise_concordance.tsv   \\
        --gtc_dir     ${params.outdir}/gtc       \\
        --out         idat_duplicate_report.tsv  \\
        --log         pairwise_concordance.log   \\
        2>&1 | tee -a concordance_run.log

    echo "========================================"  >> concordance_run.log
    echo "GENOTYPE_CONCORDANCE Completed: \$(date)" >> concordance_run.log

    cat pairwise_concordance.log >> concordance_run.log
    mv concordance_run.log pairwise_concordance.log
    """
}