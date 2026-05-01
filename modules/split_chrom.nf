process SPLIT_CHROM {
    tag "cohort"
    publishDir "${params.outdir}/qc/split_chrom", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    tuple path("chrX.bed"), path("chrX.bim"), path("chrX.fam"), emit: chrX
    tuple path("chrY.bed"), path("chrY.bim"), path("chrY.fam"), emit: chrY, optional: true
    tuple path("autosomes.bed"), path("autosomes.bim"), path("autosomes.fam"), emit: autosomes
    path "chrom_recode.log", emit: log

    script:
    def prefix = bed.baseName
    
    """
    echo "Chromosome recoding started: \$(date)" > chrom_recode.log
    
    # Check if input files exist
    if [ ! -f ${prefix}.bim ]; then
        echo "ERROR: Input files not found!" >> chrom_recode.log
        exit 1
    fi
    
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
    echo "Extracting chrX variants..." >> chrom_recode.log
    plink2 \
        --bfile ${prefix}_recode \
        --chr X \
        --make-bed \
        --allow-extra-chr \
        --out chrX 2>&1 | tee -a chrom_recode.log
    
    # Extract chrY (may not exist)
    echo "Extracting chrY variants..." >> chrom_recode.log
    plink2 \
        --bfile ${prefix}_recode \
        --chr Y \
        --make-bed \
        --allow-extra-chr \
        --out chrY 2>&1 | tee -a chrom_recode.log
    
    # Check if chrY has any variants, if not create dummy files
    if [ ! -s chrY.bim ] || [ \$(wc -l < chrY.bim 2>/dev/null || echo 0) -eq 0 ]; then
        echo "WARNING: No chrY variants found. Creating dummy chrY files." >> chrom_recode.log
        # Create dummy chrY files (copy chrX but rename to chrY - will be empty)
        echo -e "0\tY\t0\t0\tA\tC" > chrY.bim
        touch chrY.bed
        echo -e "0\t0\t0\t0\t-9" > chrY.fam
    fi
    
    # Extract autosomes (1-22)
    echo "Extracting autosome variants..." >> chrom_recode.log
    plink2 \
        --bfile ${prefix}_recode \
        --chr 1-22 \
        --make-bed \
        --allow-extra-chr \
        --out autosomes 2>&1 | tee -a chrom_recode.log
    
    # Get counts
    N_CHRX=\$(wc -l < chrX.bim 2>/dev/null || echo 0)
    N_CHRY=\$(wc -l < chrY.bim 2>/dev/null || echo 0)
    N_AUTO=\$(wc -l < autosomes.bim 2>/dev/null || echo 0)
    
    echo "Chromosome recoding completed: \$(date)" >> chrom_recode.log
    echo "========================================" >> chrom_recode.log
    echo "Summary:" >> chrom_recode.log
    echo "  chrX variants: \${N_CHRX}" >> chrom_recode.log
    echo "  chrY variants: \${N_CHRY}" >> chrom_recode.log
    echo "  Autosome variants: \${N_AUTO}" >> chrom_recode.log
    echo "========================================" >> chrom_recode.log
    
    # Warn if no chrX variants (should not happen)
    if [ \${N_CHRX} -eq 0 ]; then
        echo "ERROR: No chrX variants found! Sex checking will not work." >> chrom_recode.log
    fi
    """
}