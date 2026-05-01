process SNP_QC {
    tag "cohort"
    publishDir "${params.outdir}/qc/snp_qc", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    tuple path("cohort_snpqc.bed"),
          path("cohort_snpqc.bim"),
          path("cohort_snpqc.fam"), emit: plink
    path "snp_qc.log",              emit: log

    script:
    def prefix = bed.baseName
    
    """
    echo "========================================" > snp_qc.log
    echo "SNP QC Started: \$(date)" >> snp_qc.log
    echo "========================================" >> snp_qc.log
    
    # Initial counts
    echo "Input SNPs   : \$(wc -l < ${prefix}.bim)" >> snp_qc.log
    echo "Input samples: \$(wc -l < ${prefix}.fam)" >> snp_qc.log
    echo "" >> snp_qc.log
    
    # Step 1: Remove duplicate SNP positions
    echo "Step 1: Removing duplicate SNPs..." >> snp_qc.log
    plink \
        --bfile ${prefix} \
        --list-duplicate-vars suppress-first \
        --allow-extra-chr \
        --out dup_snps \
        --allow-no-sex

    if [ -s dup_snps.dupvar ]; then
        N_DUP=\$(wc -l < dup_snps.dupvar)
        echo "  Found \$N_DUP duplicate SNPs" >> snp_qc.log
        plink \
            --bfile ${prefix} \
            --exclude dup_snps.dupvar \
            --make-bed \
            --allow-extra-chr \
            --out no_dups \
            --allow-no-sex
    else
        echo "  No duplicate SNPs found" >> snp_qc.log
        cp ${prefix}.bed no_dups.bed
        cp ${prefix}.bim no_dups.bim
        cp ${prefix}.fam no_dups.fam
    fi

    echo "  SNPs after dedup: \$(wc -l < no_dups.bim)" >> snp_qc.log
    echo "" >> snp_qc.log

    # Step 2: Apply SNP QC filters
    echo "Step 2: Applying QC filters..." >> snp_qc.log
    echo "  --geno ${params.snp_missing} (SNP missingness)" >> snp_qc.log
    echo "  --maf  ${params.maf} (minor allele frequency)" >> snp_qc.log
    echo "  --hwe  ${params.hwe} (Hardy-Weinberg equilibrium)" >> snp_qc.log
    
    plink \
        --bfile no_dups \
        --geno ${params.snp_missing} \
        --maf ${params.maf} \
        --hwe ${params.hwe} \
        --make-bed \
        --allow-extra-chr \
        --out cohort_snpqc \
        --allow-no-sex 2>&1 | tee -a snp_qc.log

    # Step 3: Check results
    N_SNPS_QC=\$(wc -l < cohort_snpqc.bim 2>/dev/null || echo 0)
    
    echo "" >> snp_qc.log
    echo "========================================" >> snp_qc.log
    echo "QC Results:" >> snp_qc.log
    echo "  SNPs after QC: \${N_SNPS_QC}" >> snp_qc.log
    
    if [ \${N_SNPS_QC} -eq 0 ]; then
        echo "  WARNING: No SNPs passed QC filters!" >> snp_qc.log
        echo "  Attempting with more lenient filters..." >> snp_qc.log
        
        # Try more lenient filters
        plink \
            --bfile no_dups \
            --geno 0.05 \
            --maf 0.0001 \
            --hwe 1e-4 \
            --make-bed \
            --allow-extra-chr \
            --out cohort_snpqc \
            --allow-no-sex
            
        N_SNPS_QC=\$(wc -l < cohort_snpqc.bim 2>/dev/null || echo 0)
        echo "  SNPs after lenient QC: \${N_SNPS_QC}" >> snp_qc.log
        
        if [ \${N_SNPS_QC} -eq 0 ]; then
            echo "  ERROR: Still no SNPs. Using unfiltered data as fallback!" >> snp_qc.log
            # Fallback: use deduplicated but unfiltered data
            cp no_dups.bed cohort_snpqc.bed
            cp no_dups.bim cohort_snpqc.bim
            cp no_dups.fam cohort_snpqc.fam
            N_SNPS_QC=\$(wc -l < cohort_snpqc.bim)
            echo "  Using unfiltered SNPs: \${N_SNPS_QC}" >> snp_qc.log
        fi
    fi
    
    echo "========================================" >> snp_qc.log
    echo "Final summary:" >> snp_qc.log
    echo "  Samples: \$(wc -l < cohort_snpqc.fam)" >> snp_qc.log
    echo "  SNPs   : \${N_SNPS_QC}" >> snp_qc.log
    echo "========================================" >> snp_qc.log
    
    # Check if we have enough SNPs for PCA
    if [ \${N_SNPS_QC} -lt 100 ]; then
        echo "WARNING: Only \${N_SNPS_QC} SNPs available. PCA may be unreliable." >> snp_qc.log
    fi
    """
}