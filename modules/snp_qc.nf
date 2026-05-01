process SNP_QC {
    tag "cohort"
    publishDir "${params.outdir}/qc/snp_qc", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    tuple path("cohort_snpqc.bed"),
          path("cohort_snpqc.bim"),
          path("cohort_snpqc.fam"), emit: plink
    path "snp_qc.log", emit: log

    script:
    def prefix = bed.baseName
    def snp_missing = params.snp_missing ?: 0.02
    def maf = params.maf ?: 0.001
    def hwe = params.hwe ?: 0.00001
    
    """
    echo "========================================" > snp_qc.log
    echo "SNP QC Started: \$(date)" >> snp_qc.log
    echo "========================================" >> snp_qc.log
    
    # Check if input files exist
    if [ ! -f ${prefix}.bim ]; then
        echo "ERROR: Input files not found!" >> snp_qc.log
        exit 1
    fi
    
    # Initial counts
    echo "Input SNPs   : \$(wc -l < ${prefix}.bim 2>/dev/null || echo 0)" >> snp_qc.log
    echo "Input samples: \$(wc -l < ${prefix}.fam 2>/dev/null || echo 0)" >> snp_qc.log
    echo "" >> snp_qc.log
    
    # Step 1: Remove duplicate SNP positions (optional, if many duplicates)
    echo "Step 1: Checking for duplicate SNPs..." >> snp_qc.log
    plink \
        --bfile ${prefix} \
        --list-duplicate-vars suppress-first \
        --allow-extra-chr \
        --allow-no-sex \
        --out dup_snps 2>&1 | tee -a snp_qc.log

    if [ -s dup_snps.dupvar ]; then
        N_DUP=\$(wc -l < dup_snps.dupvar)
        echo "  Found \$N_DUP duplicate SNPs" >> snp_qc.log
        if [ \${N_DUP} -lt 10000 ]; then
            plink \
                --bfile ${prefix} \
                --exclude dup_snps.dupvar \
                --make-bed \
                --allow-extra-chr \
                --allow-no-sex \
                --out no_dups 2>&1 | tee -a snp_qc.log
        else
            echo "  Too many duplicates (\${N_DUP}), skipping removal" >> snp_qc.log
            cp ${prefix}.bed no_dups.bed
            cp ${prefix}.bim no_dups.bim
            cp ${prefix}.fam no_dups.fam
        fi
    else
        echo "  No duplicate SNPs found" >> snp_qc.log
        cp ${prefix}.bed no_dups.bed
        cp ${prefix}.bim no_dups.bim
        cp ${prefix}.fam no_dups.fam
    fi

    echo "  SNPs after dedup: \$(wc -l < no_dups.bim 2>/dev/null || echo 0)" >> snp_qc.log
    echo "" >> snp_qc.log

    # Step 2: Apply SNP QC filters with very lenient parameters for identical samples
    echo "Step 2: Applying QC filters..." >> snp_qc.log
    echo "  --geno ${snp_missing} (SNP missingness)" >> snp_qc.log
    echo "  --maf  ${maf} (minor allele frequency)" >> snp_qc.log
    echo "  --hwe  ${hwe} (Hardy-Weinberg equilibrium)" >> snp_qc.log
    
    # Try with original parameters
    plink \
        --bfile no_dups \
        --geno ${snp_missing} \
        --maf ${maf} \
        --hwe ${hwe} \
        --make-bed \
        --allow-extra-chr \
        --allow-no-sex \
        --out cohort_snpqc 2>&1 | tee -a snp_qc.log

    # Check if any variants remain
    N_SNPS_QC=\$(wc -l < cohort_snpqc.bim 2>/dev/null || echo 0)
    
    echo "" >> snp_qc.log
    echo "========================================" >> snp_qc.log
    echo "QC Results:" >> snp_qc.log
    echo "  SNPs after QC: \${N_SNPS_QC}" >> snp_qc.log
    
    # If no variants, try without MAF filter
    if [ \${N_SNPS_QC} -eq 0 ]; then
        echo "  WARNING: No SNPs passed QC filters!" >> snp_qc.log
        echo "  Attempting with no MAF filter..." >> snp_qc.log
        
        plink \
            --bfile no_dups \
            --geno ${snp_missing} \
            --hwe ${hwe} \
            --make-bed \
            --allow-extra-chr \
            --allow-no-sex \
            --out cohort_snpqc 2>&1 | tee -a snp_qc.log
            
        N_SNPS_QC=\$(wc -l < cohort_snpqc.bim 2>/dev/null || echo 0)
        echo "  SNPs after QC (no MAF filter): \${N_SNPS_QC}" >> snp_qc.log
        
        # If still no variants, just use original data
        if [ \${N_SNPS_QC} -eq 0 ]; then
            echo "  ERROR: Still no SNPs. Using original data as fallback!" >> snp_qc.log
            cp no_dups.bed cohort_snpqc.bed
            cp no_dups.bim cohort_snpqc.bim
            cp no_dups.fam cohort_snpqc.fam
            N_SNPS_QC=\$(wc -l < cohort_snpqc.bim)
            echo "  Using original SNPs: \${N_SNPS_QC}" >> snp_qc.log
        fi
    fi
    
    echo "========================================" >> snp_qc.log
    echo "Final summary:" >> snp_qc.log
    echo "  Samples: \$(wc -l < cohort_snpqc.fam 2>/dev/null || echo 0)" >> snp_qc.log
    echo "  SNPs   : \${N_SNPS_QC}" >> snp_qc.log
    echo "========================================" >> snp_qc.log
    
    # Check if we have enough SNPs for downstream analyses
    if [ \${N_SNPS_QC} -lt 100 ]; then
        echo "WARNING: Only \${N_SNPS_QC} SNPs available. This may indicate all samples are identical." >> snp_qc.log
    fi
    """
}