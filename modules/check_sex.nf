process CHECK_SEX {
    publishDir "${params.outdir}/qc/sex_check", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)
    path sex_info

    output:
    path "sex_check.sexcheck",          emit: sexcheck
    path "sex_check_annotated.tsv",     emit: annotated
    path "sex_discordant_samples.txt",  emit: discordant

    script:
    """
    PREFIX=\$(basename ${bed} .bed)

    # ── PLINK sex check using chrX F-statistic ───────────────────────────
    # Females: F < 0.2   Males: F > 0.8   (PLINK defaults)
    plink \
        --bfile \${PREFIX} \
        --check-sex \
        --allow-extra-chr \
        --out sex_check

    # ── Annotate with collected gender and flag discordant samples ────────
    python ${projectDir}/bin/annotate_sex_check.py \
        --sexcheck   sex_check.sexcheck \
        --sex_info   ${sex_info}         \
        --out_annot  sex_check_annotated.tsv \
        --out_discord sex_discordant_samples.txt
    """
}
