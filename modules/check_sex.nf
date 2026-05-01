process CHECK_SEX {
    tag "cohort"
    publishDir "${params.outdir}/qc/check_sex", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)  // chrX from SPLIT_CHROM
    path sex_info  // from PREP_INPUTS (collected sex)

    output:
    path "sexcheck.txt", emit: sexcheck
    path "sex_discordance.txt", emit: discordance
    path "sex_plot.png", emit: plot

    script:
    """
    echo "Sex check started: \$(date)" > sex_check.log
    
    # Step 1: Run PLINK's sex check on chrX
    # This calculates F-statistic (inbreeding coefficient on X chromosome)
    # F < 0.2 = female, F > 0.8 = male, 0.2-0.8 = ambiguous
    plink2 \
        --bfile ${bed.baseName} \
        --check-sex \
        --allow-extra-chr \
        --out sexcheck
    
    # Step 2: Parse results and compare with collected sex
    python ${projectDir}/bin/compare_sex.py \
        --plink_sex sexcheck.sexcheck \
        --collected_sex ${sex_info} \
        --out sex_discordance.txt \
        --plot sex_plot.png
    
    echo "Sex check completed: \$(date)" >> sex_check.log
    
    # Summary statistics
    echo "" >> sex_check.log
    echo "Summary:" >> sex_check.log
    echo "  Total samples: \$(wc -l < sexcheck.sexcheck)" >> sex_check.log
    echo "  Discordant samples: \$(grep -c PROBLEM sex_discordance.txt 2>/dev/null || echo 0)" >> sex_check.log
    """
}