process PCA {
    publishDir "${params.outdir}/qc/pca", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)
    path keep_list

    output:
    path "pca.eigenvec",    emit: eigenvec
    path "pca.eigenval",    emit: eigenval
    path "pca_plot.png",    emit: plot

    script:
    """
    PREFIX=\$(basename ${bed} .bed)

    # ── 1. Keep QC-passing samples + LD pruning ───────────────────────────
    plink2 \
        --bfile \${PREFIX} \
        --keep ${keep_list} \
        --maf ${params.maf} \
        --indep-pairwise 200 50 0.25 \
        --allow-extra-chr \
        --out pruned

    plink2 \
        --bfile \${PREFIX} \
        --keep ${keep_list} \
        --extract pruned.prune.in \
        --make-bed \
        --allow-extra-chr \
        --out pruned_data

    # ── 2. PCA ────────────────────────────────────────────────────────────
    plink2 \
        --bfile pruned_data \
        --pca ${params.pca_components} \
        --allow-extra-chr \
        --out pca

    # ── 3. Plot ───────────────────────────────────────────────────────────
    python ${projectDir}/bin/plot_pca.py \
        --eigenvec pca.eigenvec \
        --eigenval pca.eigenval \
        --out      pca_plot.png
    """
}
