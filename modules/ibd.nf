process IBD {
    publishDir "${params.outdir}/qc/ibd", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)
    path keep_list

    output:
    path "ibd.genome",              emit: genome
    path "ibd_duplicates.tsv",      emit: duplicates
    path "ibd_plot.png",            emit: plot

    script:
    """
    PREFIX=\$(basename ${bed} .bed)

    # ── 1. Keep only QC-passing samples ───────────────────────────────────
    plink \
        --bfile \${PREFIX} \
        --keep ${keep_list} \
        --make-bed \
        --allow-extra-chr \
        --out qc_passed

    # ── 2. LD pruning (recommended before IBD) ────────────────────────────
    plink \
        --bfile qc_passed \
        --maf ${params.maf} \
        --indep-pairwise 50 10 0.1 \
        --allow-extra-chr \
        --out pruned

    plink \
        --bfile qc_passed \
        --extract pruned.prune.in \
        --make-bed \
        --allow-extra-chr \
        --out pruned_data

    # ── 3. IBD / IBS estimation ───────────────────────────────────────────
    plink \
        --bfile pruned_data \
        --genome \
        --min ${params.ibd_pi_hat} \
        --allow-extra-chr \
        --out ibd

    # ── 4. Flag duplicates and cross-reference with sex discordance ───────
    python ${projectDir}/bin/flag_ibd_duplicates.py \
        --genome   ibd.genome \
        --pi_hat   ${params.ibd_pi_hat} \
        --out      ibd_duplicates.tsv  \
        --plot     ibd_plot.png
    """
}
