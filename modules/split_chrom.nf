process SPLIT_CHROM {
    publishDir "${params.outdir}/qc/chromosomes", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    tuple path("chrX.bed"), path("chrX.bim"), path("chrX.fam"), emit: chrX
    tuple path("chrY.bed"), path("chrY.bim"), path("chrY.fam"), emit: chrY

    script:
    """
    PREFIX=\$(basename ${bed} .bed)

    # Extract chrX (handles both '23' and 'X' chromosome encoding)
    plink \
        --bfile \${PREFIX} \
        --chr 23 \
        --make-bed \
        --allow-extra-chr \
        --out chrX

    # Extract chrY
    plink \
        --bfile \${PREFIX} \
        --chr 24 \
        --make-bed \
        --allow-extra-chr \
        --out chrY
    """
}
