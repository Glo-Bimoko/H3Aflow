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
    def prefix      = bed.baseName
    def snp_missing = params.snp_missing ?: 0.05
    def maf         = params.maf         ?: 0.0001
    def hwe         = params.hwe         ?: 1e-4

    """
    #!/usr/bin/env bash
    set -euo pipefail

    # ════════════════════════════════════════════════════════════════════════
    # SNP QC — autosomes only
    #
    # chrX, chrY, and XY (PAR) are excluded from --geno / --maf / --hwe.
    # Sex-stratified chrX QC is handled entirely inside CHECK_SEX, which
    # runs in parallel from the same FILTER_SAMPLES output.  Applying
    # cohort-wide filters to chrX before sex check is invalid because:
    #
    #   --geno : males are hemizygous on non-PAR X; het-haploid calls inflate
    #            per-SNP missingness, removing the most informative female-
    #            heterozygous X SNPs.
    #   --hwe  : HWE is undefined for hemizygous loci; spurious p-values on a
    #            mixed-sex cohort drop good non-PAR X SNPs.
    #   --maf  : cohort-wide MAF conflates male hemizygous and female diploid
    #            allele counts, making X-specific rare variants appear absent.
    #
    # The previous version split autosomes and sex chromosomes into separate
    # PLINK datasets and merged them back with --bmerge.  bmerge was silently
    # dropping autosomal SNPs due to chromosome code handling.  The --not-chr
    # flag avoids the split/merge entirely and is safer.
    # ════════════════════════════════════════════════════════════════════════

    LOG=snp_qc.log
    echo "========================================"  > \$LOG
    echo "SNP QC Started: \$(date)"                 >> \$LOG
    echo "========================================"  >> \$LOG
    echo ""                                          >> \$LOG
    echo "Parameters:"                               >> \$LOG
    echo "  --geno (snp_missing) : ${snp_missing}"  >> \$LOG
    echo "  --maf                : ${maf}"           >> \$LOG
    echo "  --hwe                : ${hwe}"           >> \$LOG
    echo "  Chromosomes excluded : X, Y, XY (PAR)"  >> \$LOG
    echo ""                                          >> \$LOG

    if [ ! -f ${prefix}.bim ]; then
        echo "ERROR: Input BIM not found!" >> \$LOG
        exit 1
    fi

    N_IN=\$(wc -l < ${prefix}.bim)
    N_SAMPLES=\$(wc -l < ${prefix}.fam)
    echo "Input variants : \${N_IN}"   >> \$LOG
    echo "Input samples  : \${N_SAMPLES}" >> \$LOG
    echo "Chromosome codes in input:"   >> \$LOG
    awk '{print \$1}' ${prefix}.bim | sort -u >> \$LOG
    echo "" >> \$LOG

    # ── Step 1: Remove duplicate SNP positions ────────────────────────────────
    # Kept from the previous snp_qc.nf — duplicate positions cause PLINK2
    # and downstream tools to fail; removing them here is valid regardless
    # of chromosome and unrelated to the sex-chromosome filtering change.
    echo "--- Step 1: Remove duplicate SNPs ---" >> \$LOG
    plink \\
        --bfile ${prefix} \\
        --list-duplicate-vars suppress-first \\
        --allow-extra-chr \\
        --allow-no-sex \\
        --out dup_snps >> \$LOG 2>&1 || true

    if [ -s dup_snps.dupvar ]; then
        N_DUP=\$(wc -l < dup_snps.dupvar)
        echo "  Removing \${N_DUP} duplicate SNPs..." >> \$LOG
        plink \\
            --bfile ${prefix} \\
            --exclude dup_snps.dupvar \\
            --make-bed \\
            --allow-extra-chr \\
            --allow-no-sex \\
            --out no_dups >> \$LOG 2>&1
    else
        echo "  No duplicates found." >> \$LOG
        ln -sf ${prefix}.bed no_dups.bed
        ln -sf ${prefix}.bim no_dups.bim
        ln -sf ${prefix}.fam no_dups.fam
    fi
    echo "  Variants after dedup: \$(wc -l < no_dups.bim)" >> \$LOG
    echo "" >> \$LOG

    # ── Step 2: Autosomal SNP QC (--not-chr X Y XY) ──────────────────────────
    # --not-chr tells PLINK to apply all subsequent filters only to the
    # specified chromosomes' complement.  chrX/Y/XY SNPs are silently excluded
    # from the output — no split, no merge, no bmerge corruption risk.
    echo "--- Step 2: Autosomal SNP QC (--not-chr X Y XY) ---" >> \$LOG

    N_X=\$(awk '\$1==23 || \$1=="X"'   no_dups.bim | wc -l || echo 0)
    N_Y=\$(awk '\$1==24 || \$1=="Y"'   no_dups.bim | wc -l || echo 0)
    N_XY=\$(awk '\$1==25 || \$1=="XY"' no_dups.bim | wc -l || echo 0)
    N_AUTO_IN=\$(( \$(wc -l < no_dups.bim) - N_X - N_Y - N_XY ))
    echo "  Autosomal SNPs (input)    : \${N_AUTO_IN}" >> \$LOG
    echo "  chrX SNPs (excluded)      : \${N_X}"       >> \$LOG
    echo "  chrY SNPs (excluded)      : \${N_Y}"       >> \$LOG
    echo "  chrXY/PAR SNPs (excluded) : \${N_XY}"      >> \$LOG
    echo "" >> \$LOG

    plink \\
        --bfile no_dups \\
        --not-chr X Y XY \\
        --geno ${snp_missing} \\
        --maf  ${maf} \\
        --hwe  ${hwe} \\
        --make-bed \\
        --allow-no-sex \\
        --out cohort_snpqc >> \$LOG 2>&1

    N_OUT=\$(wc -l < cohort_snpqc.bim 2>/dev/null || echo 0)

    # Safety net: if MAF removes everything (can happen on tiny pilot datasets),
    # retry without --maf and warn loudly.
    if [ "\${N_OUT}" -eq 0 ]; then
        echo "  WARNING: All autosomal SNPs removed — retrying without --maf..." >> \$LOG
        plink \\
            --bfile no_dups \\
            --not-chr X Y XY \\
            --geno ${snp_missing} \\
            --hwe  ${hwe} \\
            --make-bed \\
            --allow-no-sex \\
            --out cohort_snpqc >> \$LOG 2>&1
        N_OUT=\$(wc -l < cohort_snpqc.bim 2>/dev/null || echo 0)
        echo "  Autosomal SNPs (no MAF fallback): \${N_OUT}" >> \$LOG
    fi

    N_DROPPED=\$(( N_AUTO_IN - N_OUT ))

    # ── Final summary ─────────────────────────────────────────────────────────
    echo "" >> \$LOG
    echo "Chromosome codes in final dataset:" >> \$LOG
    awk '{print \$1}' cohort_snpqc.bim | sort -u >> \$LOG
    echo "" >> \$LOG
    echo "========================================"     >> \$LOG
    echo "SNP QC Summary"                               >> \$LOG
    echo "========================================"     >> \$LOG
    echo "  Input variants (all chr)    : \${N_IN}"     >> \$LOG
    echo "  After dedup                 : \$(wc -l < no_dups.bim)" >> \$LOG
    echo "  Autosomal input             : \${N_AUTO_IN}" >> \$LOG
    echo "  Autosomal after QC          : \${N_OUT}"    >> \$LOG
    echo "  Autosomal dropped           : \${N_DROPPED}" >> \$LOG
    echo "  chrX/Y/XY excluded from QC  : \$(( N_X + N_Y + N_XY )) (sex check uses GTC computed_gender)" >> \$LOG
    echo "  Output samples              : \$(wc -l < cohort_snpqc.fam)" >> \$LOG
    echo "========================================"     >> \$LOG
    echo "SNP QC Completed: \$(date)"                   >> \$LOG

    if [ "\${N_OUT}" -eq 0 ]; then
        echo "ERROR: Zero autosomal SNPs in final dataset — cannot continue." >> \$LOG
        exit 1
    fi
    """
}