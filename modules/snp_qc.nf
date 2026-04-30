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
    # Step 1: Remove duplicate SNP positions
    plink \
        --bfile  ${prefix}          \
        --list-duplicate-vars suppress-first \
        --out    dup_snps           \
        --allow-extra-chr

    if [ -s dup_snps.dupvar ]; then
        plink \
            --bfile        ${prefix}    \
            --exclude      dup_snps.dupvar \
            --make-bed     \
            --out          no_dups      \
            --allow-extra-chr
    else
        cp ${prefix}.bed no_dups.bed
        cp ${prefix}.bim no_dups.bim
        cp ${prefix}.fam no_dups.fam
    fi

    # Step 2: SNP-level QC
    #   --geno ${params.snp_missing}  remove SNPs with missingness > threshold
    #   --maf  ${params.maf}          remove SNPs with MAF below threshold
    #   --hwe  ${params.hwe}          remove SNPs failing HWE p-value threshold
    plink \
        --bfile        no_dups         \
        --geno         ${params.snp_missing} \
        --maf          ${params.maf}   \
        --hwe          ${params.hwe}   \
        --make-bed     \
        --out          cohort_snpqc    \
        --allow-extra-chr             \
        --allow-no-sex

    # Collect stats
    echo "SNP QC summary" > snp_qc.log
    echo "Input SNPs   : \$(wc -l < ${prefix}.bim)"     >> snp_qc.log
    echo "Duplicate SNPs removed: \$(wc -l < dup_snps.dupvar 2>/dev/null || echo 0)" >> snp_qc.log
    echo "Output SNPs  : \$(wc -l < cohort_snpqc.bim)"  >> snp_qc.log
    echo "Output samples: \$(wc -l < cohort_snpqc.fam)" >> snp_qc.log
    """
}
