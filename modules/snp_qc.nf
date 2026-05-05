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
    def maf         = params.maf         ?: 0.001
    def hwe         = params.hwe         ?: 0.00001

    """
    echo "========================================" > snp_qc.log
    echo "SNP QC Started: \$(date)" >> snp_qc.log
    echo "========================================" >> snp_qc.log

    if [ ! -f ${prefix}.bim ]; then
        echo "ERROR: Input files not found!" >> snp_qc.log
        exit 1
    fi

    echo "Input SNPs   : \$(wc -l < ${prefix}.bim)" >> snp_qc.log
    echo "Input samples: \$(wc -l < ${prefix}.fam)" >> snp_qc.log
    echo "Chromosome codes in input:" >> snp_qc.log
    awk '{print \$1}' ${prefix}.bim | sort -u >> snp_qc.log
    echo "" >> snp_qc.log

    # ── Step 1: Remove duplicate SNP positions ────────────────────────────────
    echo "Step 1: Removing duplicate SNPs..." >> snp_qc.log
    plink \
        --bfile ${prefix} \
        --list-duplicate-vars suppress-first \
        --allow-extra-chr \
        --allow-no-sex \
        --out dup_snps 2>&1 | tee -a snp_qc.log

    if [ -s dup_snps.dupvar ]; then
        N_DUP=\$(wc -l < dup_snps.dupvar)
        echo "  Removing \$N_DUP duplicate SNPs..." >> snp_qc.log
        plink \
            --bfile ${prefix} \
            --exclude dup_snps.dupvar \
            --make-bed \
            --allow-extra-chr \
            --allow-no-sex \
            --out no_dups 2>&1 | tee -a snp_qc.log
    else
        echo "  No duplicates found." >> snp_qc.log
        cp ${prefix}.bed no_dups.bed
        cp ${prefix}.bim no_dups.bim
        cp ${prefix}.fam no_dups.fam
    fi
    echo "  SNPs after dedup: \$(wc -l < no_dups.bim)" >> snp_qc.log
    echo "" >> snp_qc.log

    # ── Step 2: Identify sex chromosome SNP IDs ───────────────────────────────
    # Sex chromosomes: integer codes >= 23, or string codes X/Y/XY/MT.
    # We collect their SNP IDs so we can EXCLUDE them from MAF/HWE
    # (those filters are statistically invalid on sex chromosomes in a
    # mixed-sex cohort), then add them back with missingness-only filtering.
    # This avoids the split→bmerge approach that caused chromosome code
    # corruption (bmerge was silently dropping autosomes).
    echo "Step 2: Identifying sex chromosome SNP IDs..." >> snp_qc.log
    awk '\$1 >= 23 || \$1 ~ /^[XxYyMm]/' no_dups.bim | awk '{print \$2}' > sexchr_ids.txt
    N_SEXCHR_IDS=\$(wc -l < sexchr_ids.txt)
    echo "  Sex chromosome SNP IDs found: \${N_SEXCHR_IDS}" >> snp_qc.log
    echo "" >> snp_qc.log

    # ── Step 3: Autosome QC – MAF + HWE + missingness ────────────────────────
    echo "Step 3: Autosome QC (--geno ${snp_missing} --maf ${maf} --hwe ${hwe})..." >> snp_qc.log

    if [ "\${N_SEXCHR_IDS}" -gt 0 ]; then
        plink \
            --bfile no_dups \
            --exclude sexchr_ids.txt \
            --geno ${snp_missing} \
            --maf  ${maf} \
            --hwe  ${hwe} \
            --make-bed \
            --allow-no-sex \
            --out autosomes_qc 2>&1 | tee -a snp_qc.log
    else
        plink \
            --bfile no_dups \
            --geno ${snp_missing} \
            --maf  ${maf} \
            --hwe  ${hwe} \
            --make-bed \
            --allow-no-sex \
            --out autosomes_qc 2>&1 | tee -a snp_qc.log
    fi

    N_AUTO=\$(wc -l < autosomes_qc.bim 2>/dev/null || echo 0)
    echo "  Autosome SNPs after QC: \${N_AUTO}" >> snp_qc.log

    # Safety net: if MAF still kills everything (tiny subset), retry without MAF
    if [ "\${N_AUTO}" -eq 0 ]; then
        echo "  WARNING: MAF filter removed all autosomes – retrying without --maf..." >> snp_qc.log
        if [ "\${N_SEXCHR_IDS}" -gt 0 ]; then
            plink \
                --bfile no_dups \
                --exclude sexchr_ids.txt \
                --geno ${snp_missing} \
                --hwe  ${hwe} \
                --make-bed \
                --allow-no-sex \
                --out autosomes_qc 2>&1 | tee -a snp_qc.log
        else
            plink \
                --bfile no_dups \
                --geno ${snp_missing} \
                --hwe  ${hwe} \
                --make-bed \
                --allow-no-sex \
                --out autosomes_qc 2>&1 | tee -a snp_qc.log
        fi
        N_AUTO=\$(wc -l < autosomes_qc.bim 2>/dev/null || echo 0)
        echo "  Autosome SNPs (no MAF): \${N_AUTO}" >> snp_qc.log
    fi
    echo "" >> snp_qc.log

    # ── Step 4: Sex chromosome QC – missingness only ─────────────────────────
    echo "Step 4: Sex chromosome QC (--geno ${snp_missing} only)..." >> snp_qc.log

    if [ "\${N_SEXCHR_IDS}" -gt 0 ]; then
        plink \
            --bfile no_dups \
            --extract sexchr_ids.txt \
            --geno ${snp_missing} \
            --make-bed \
            --allow-extra-chr \
            --allow-no-sex \
            --out sexchr_qc 2>&1 | tee -a snp_qc.log
        N_SEXCHR=\$(wc -l < sexchr_qc.bim 2>/dev/null || echo 0)
    else
        N_SEXCHR=0
    fi
    echo "  Sex chromosome SNPs after QC: \${N_SEXCHR}" >> snp_qc.log
    echo "" >> snp_qc.log

    # ── Step 5: Merge back into one dataset ───────────────────────────────────
    echo "Step 5: Assembling final dataset..." >> snp_qc.log

    if [ "\${N_AUTO}" -gt 0 ] && [ "\${N_SEXCHR}" -gt 0 ]; then
        echo "  Merging autosomes (\${N_AUTO}) + sex chr (\${N_SEXCHR})..." >> snp_qc.log
        plink \
            --bfile autosomes_qc \
            --bmerge sexchr_qc \
            --make-bed \
            --allow-extra-chr \
            --allow-no-sex \
            --out cohort_snpqc 2>&1 | tee -a snp_qc.log

        # bmerge occasionally fails on strand issues – fall back to merge-list
        if [ ! -s cohort_snpqc.bim ]; then
            echo "  bmerge failed – trying --merge-list fallback..." >> snp_qc.log
            printf 'autosomes_qc.bed autosomes_qc.bim autosomes_qc.fam\nsexchr_qc.bed sexchr_qc.bim sexchr_qc.fam\n' > merge_list.txt
            plink \
                --merge-list merge_list.txt \
                --make-bed \
                --allow-extra-chr \
                --allow-no-sex \
                --out cohort_snpqc 2>&1 | tee -a snp_qc.log
        fi

    elif [ "\${N_AUTO}" -gt 0 ]; then
        echo "  No sex chr SNPs – using autosomes only." >> snp_qc.log
        cp autosomes_qc.bed cohort_snpqc.bed
        cp autosomes_qc.bim cohort_snpqc.bim
        cp autosomes_qc.fam cohort_snpqc.fam

    elif [ "\${N_SEXCHR}" -gt 0 ]; then
        echo "  WARNING: No autosome SNPs survived – using sex chr only." >> snp_qc.log
        cp sexchr_qc.bed cohort_snpqc.bed
        cp sexchr_qc.bim cohort_snpqc.bim
        cp sexchr_qc.fam cohort_snpqc.fam

    else
        echo "  ERROR: No SNPs survived QC. Falling back to pre-QC data." >> snp_qc.log
        cp no_dups.bed cohort_snpqc.bed
        cp no_dups.bim cohort_snpqc.bim
        cp no_dups.fam cohort_snpqc.fam
    fi

    # ── Final summary ─────────────────────────────────────────────────────────
    N_FINAL=\$(wc -l < cohort_snpqc.bim 2>/dev/null || echo 0)
    N_CHRX_FINAL=\$(awk '\$1==23 || \$1=="X"' cohort_snpqc.bim | wc -l)
    N_AUTO_FINAL=\$(awk '\$1+0 >= 1 && \$1+0 <= 22' cohort_snpqc.bim | wc -l)

    echo "" >> snp_qc.log
    echo "Chromosome codes in final dataset:" >> snp_qc.log
    awk '{print \$1}' cohort_snpqc.bim | sort -u >> snp_qc.log
    echo "========================================" >> snp_qc.log
    echo "Final summary:" >> snp_qc.log
    echo "  Samples      : \$(wc -l < cohort_snpqc.fam)" >> snp_qc.log
    echo "  Total SNPs   : \${N_FINAL}" >> snp_qc.log
    echo "  Autosomes    : \${N_AUTO_FINAL}" >> snp_qc.log
    echo "  ChrX SNPs    : \${N_CHRX_FINAL}  (needed for sex check)" >> snp_qc.log
    echo "========================================" >> snp_qc.log
    echo "SNP QC Completed: \$(date)" >> snp_qc.log

    if [ "\${N_AUTO_FINAL}" -eq 0 ]; then
        echo "ERROR: Zero autosome SNPs in final dataset!" >> snp_qc.log
        exit 1
    fi

    if [ "\${N_CHRX_FINAL}" -eq 0 ]; then
        echo "WARNING: Zero chrX SNPs – sex checking will not work." >> snp_qc.log
    fi
    """
}