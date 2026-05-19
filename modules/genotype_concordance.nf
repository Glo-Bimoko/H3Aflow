process GENOTYPE_CONCORDANCE {
    tag "cohort"
    publishDir "${params.outdir}/qc/concordance", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    path "pairwise_concordance.tsv", emit: concordance
    path "pairwise_concordance.log", emit: log

    script:
    def prefix        = bed.baseName
    def chunk         = params.concordance_chunk     ?: 50000
    def dup_threshold = params.concordance_dup_thresh ?: 99.5

    """
    echo "========================================" > concordance_run.log
    echo "GENOTYPE_CONCORDANCE Started: \$(date)"  >> concordance_run.log
    echo "========================================"  >> concordance_run.log
    echo "Input prefix : ${prefix}"                 >> concordance_run.log
    echo "Samples      : \$(wc -l < ${fam})"        >> concordance_run.log
    echo "SNPs         : \$(wc -l < ${bim})"        >> concordance_run.log

    # ── Step 1: recode BED → .traw (SNP-major, A-transpose) ──────────────────
    # .traw format: one ROW per SNP, one COLUMN per sample.
    # This lets pairwise_concordance.py stream SNP chunks without ever
    # loading the full matrix — peak memory stays flat regardless of N.
    #
    # --recode A-transpose : additive encoding (0/1/2); missing = NA
    # --autosome-num 26    : treat codes up to 26 as autosomal (handles
    #                        chrXY/PAR coded as 25/26 in this dataset)
    # --allow-extra-chr    : pass non-standard contigs through silently
    # --allow-no-sex       : skip the sex-check warning (sex is set later)
    plink \\
        --bfile   ${prefix}        \\
        --recode  A-transpose      \\
        --autosome-num 26          \\
        --allow-extra-chr          \\
        --allow-no-sex             \\
        --out     cohort_for_conc  \\
        2>&1 | tee -a concordance_run.log

    echo "  .traw written: cohort_for_conc.traw" >> concordance_run.log

    # ── Step 2: pairwise concordance (chunked streaming) ─────────────────────
    python3 ${projectDir}/bin/pairwise_concordance.py \\
        --traw          cohort_for_conc.traw     \\
        --out           pairwise_concordance.tsv \\
        --chunk         ${chunk}                 \\
        --dup-threshold ${dup_threshold}          \\
        2>&1 | tee -a concordance_run.log

    echo "========================================"  >> concordance_run.log
    echo "GENOTYPE_CONCORDANCE Completed: \$(date)" >> concordance_run.log

    # Append the Python-generated summary to the main log
    cat pairwise_concordance.log >> concordance_run.log
    mv concordance_run.log pairwise_concordance.log
    """
}
