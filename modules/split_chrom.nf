process SPLIT_CHROM {
    tag "cohort"
    publishDir "${params.outdir}/qc/split_chrom", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    tuple path("chrX.bed"), path("chrX.bim"), path("chrX.fam"), emit: chrX
    tuple path("chrY.bed"), path("chrY.bim"), path("chrY.fam"), emit: chrY
    tuple path("autosomes.bed"), path("autosomes.bim"), path("autosomes.fam"), emit: autosomes
    path "chrom_recode.log", emit: log

    script:
    def prefix = bed.baseName
    
    """
    echo "Chromosome recoding started: \$(date)" > chrom_recode.log
    
    # Create a recoded version with standard chromosome names
    # Map 23->X, 24->Y, 25->XY/MT (if needed)
    awk '{
        if (\$1 == 23) \$1 = "X"
        else if (\$1 == 24) \$1 = "Y"
        else if (\$1 == 25) \$1 = "XY"
        print \$0
    }' ${prefix}.bim > ${prefix}_recode.bim
    
    # Keep original files but with recoded bim
    cp ${prefix}.bed ${prefix}_recode.bed
    cp ${prefix}.fam ${prefix}_recode.fam
    
    # Extract chrX (now encoded as "X")
    plink2 \
        --bfile ${prefix}_recode \
        --chr X \
        --make-bed \
        --allow-extra-chr \
        --out chrX
    
    # Extract chrY
    plink2 \
        --bfile ${prefix}_recode \
        --chr Y \
        --make-bed \
        --allow-extra-chr \
        --out chrY
    
    # Extract autosomes (1-22)
    plink2 \
        --bfile ${prefix}_recode \
        --chr 1-22 \
        --make-bed \
        --allow-extra-chr \
        --out autosomes
    
    echo "Chromosome recoding completed: \$(date)" >> chrom_recode.log
    echo "  chrX variants: \$(wc -l < chrX.bim 2>/dev/null || echo 0)" >> chrom_recode.log
    echo "  chrY variants: \$(wc -l < chrY.bim 2>/dev/null || echo 0)" >> chrom_recode.log
    echo "  Autosome variants: \$(wc -l < autosomes.bim 2>/dev/null || echo 0)" >> chrom_recode.log
    """
}