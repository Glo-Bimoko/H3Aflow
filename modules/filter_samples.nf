process FILTER_SAMPLES {
    tag "cohort"
    publishDir "${params.outdir}/qc/sample_qc", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)
    path keep_list

    output:
    tuple path("cohort_filtered.bed"),
          path("cohort_filtered.bim"),
          path("cohort_filtered.fam"), emit: plink

    script:
    def prefix = bed.baseName
    """
    n_keep=\$(wc -l < ${keep_list})
    echo "[filter_samples] Keeping \${n_keep} samples from ${prefix}"

    plink \
        --bfile        ${prefix}    \
        --keep         ${keep_list} \
        --make-bed                  \
        --out          cohort_filtered \
        --allow-extra-chr           \
        --allow-no-sex
    """
}
