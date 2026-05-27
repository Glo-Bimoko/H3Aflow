process SPLIT_CHROM {
    tag "cohort"
    publishDir "${params.outdir}/qc/split_chrom", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    tuple path("chrX.bed"),      path("chrX.bim"),      path("chrX.fam"),      emit: chrX
    tuple path("chrY.bed"),      path("chrY.bim"),      path("chrY.fam"),      emit: chrY,      optional: true
    tuple path("autosomes.bed"), path("autosomes.bim"), path("autosomes.fam"), emit: autosomes
    path "chrom_recode.log", emit: log

    script:
    def prefix = bed.baseName

    """
    echo "========================================" > chrom_recode.log
    echo "SPLIT_CHROM Started: \$(date)" >> chrom_recode.log
    echo "========================================" >> chrom_recode.log

    if [ ! -f ${prefix}.bim ]; then
        echo "ERROR: Input bim not found!" >> chrom_recode.log
        exit 1
    fi

    TOTAL_IN=\$(wc -l < ${prefix}.bim)
    echo "Input variants : \${TOTAL_IN}" >> chrom_recode.log
    echo "Unique chr codes in input:" >> chrom_recode.log
    awk '{print \$1}' ${prefix}.bim | sort -u >> chrom_recode.log
    echo "" >> chrom_recode.log

    # ── Normalise chromosome codes in .bim ───────────────────────────────────
    # After plink --bmerge the sex-chr partition (written by snp_qc.nf with
    # --allow-extra-chr) may already have "X"/"Y" string codes while autosomes
    # are integers 1-22.  We normalise everything so downstream plink calls
    # see consistent codes:
    #   23 or X  → X
    #   24 or Y  → Y
    #   25       → XY
    #   1-22     → unchanged integers
    awk 'BEGIN{OFS="\\t"} {
        chr = \$1
        if (chr == "23" || chr == "X")  chr = "X"
        else if (chr == "24" || chr == "Y")  chr = "Y"
        else if (chr == "25" || chr == "XY") chr = "XY"
        \$1 = chr
        print
    }' ${prefix}.bim > recoded.bim

    cp ${prefix}.bed recoded.bed
    cp ${prefix}.fam recoded.fam

    echo "Unique chr codes after recode:" >> chrom_recode.log
    awk '{print \$1}' recoded.bim | sort -u >> chrom_recode.log
    echo "" >> chrom_recode.log

    # ── Extract chrX ─────────────────────────────────────────────────────────
    # Use plink 1.9: it handles mixed integer/string chr codes more gracefully
    # than plink2 when --allow-extra-chr is set.
    echo "Extracting chrX..." >> chrom_recode.log
    plink \\
        --bfile recoded \\
        --chr X \\
        --make-bed \\
        --allow-extra-chr \\
        --allow-no-sex \\
        --out chrX 2>&1 | tee -a chrom_recode.log

    N_CHRX=\$(wc -l < chrX.bim 2>/dev/null || echo 0)
    echo "  chrX variants: \${N_CHRX}" >> chrom_recode.log

    if [ "\${N_CHRX}" -eq 0 ]; then
        echo "ERROR: No chrX variants found after recode!" >> chrom_recode.log
        echo "  Check that snp_qc.nf sexchr_qc step is producing X-coded SNPs." >> chrom_recode.log
    fi

    # ── Extract chrY ─────────────────────────────────────────────────────────
    echo "Extracting chrY..." >> chrom_recode.log
    plink \\
        --bfile recoded \\
        --chr Y \\
        --make-bed \\
        --allow-extra-chr \\
        --allow-no-sex \\
        --out chrY 2>&1 | tee -a chrom_recode.log

    N_CHRY=\$(wc -l < chrY.bim 2>/dev/null || echo 0)
    echo "  chrY variants: \${N_CHRY}" >> chrom_recode.log

    # Create dummy chrY files if none found (common – array has few Y probes)
    if [ "\${N_CHRY}" -eq 0 ]; then
        echo "  WARNING: No chrY variants – creating dummy files." >> chrom_recode.log
        printf '0\tY\t0\t1\tA\tC\n' > chrY.bim
        # Minimal valid .bed: magic bytes + 1 SNP x 0 samples = just the header
        printf '\\x6c\\x1b\\x01' > chrY.bed
        # .fam: copy from input so sample list is consistent
        cp ${prefix}.fam chrY.fam
    fi

    # ── Extract autosomes (chr 1-22) ─────────────────────────────────────────
    echo "Extracting autosomes (1-22)..." >> chrom_recode.log
    plink \\
        --bfile recoded \\
        --chr 1-22 \\
        --make-bed \\
        --allow-extra-chr \\
        --allow-no-sex \\
        --out autosomes 2>&1 | tee -a chrom_recode.log

    N_AUTO=\$(wc -l < autosomes.bim 2>/dev/null || echo 0)
    echo "  Autosome variants: \${N_AUTO}" >> chrom_recode.log

    # ── Fallback: if plink 1-22 range fails, extract by excluding sex chrs ───
    # This handles edge cases where chr codes are non-standard after merge.
    if [ "\${N_AUTO}" -eq 0 ]; then
        echo "  WARNING: --chr 1-22 returned 0 variants. Trying exclusion fallback..." >> chrom_recode.log

        # Extract SNP IDs that are NOT on X, Y, XY, MT
        awk '\$1 !~ /^[XxYyMm]/ && \$1+0 >= 1 && \$1+0 <= 22 {print \$2}' recoded.bim > autosome_snps.txt
        N_AUTO_IDS=\$(wc -l < autosome_snps.txt)
        echo "  Autosome SNP IDs found by awk: \${N_AUTO_IDS}" >> chrom_recode.log

        if [ "\${N_AUTO_IDS}" -gt 0 ]; then
            plink \\
                --bfile recoded \\
                --extract autosome_snps.txt \\
                --make-bed \\
                --allow-extra-chr \\
                --allow-no-sex \\
                --out autosomes 2>&1 | tee -a chrom_recode.log
            N_AUTO=\$(wc -l < autosomes.bim 2>/dev/null || echo 0)
            echo "  Autosome variants after fallback: \${N_AUTO}" >> chrom_recode.log
        else
            echo "  ERROR: No autosome SNPs identifiable – check chromosome coding." >> chrom_recode.log
            # Create minimal dummy autosomes so Nextflow output check passes,
            # but exit 1 so the user sees the failure clearly.
            exit 1
        fi
    fi

    # ── Final summary ─────────────────────────────────────────────────────────
    echo "" >> chrom_recode.log
    echo "========================================" >> chrom_recode.log
    echo "Summary:" >> chrom_recode.log
    echo "  Input variants  : \${TOTAL_IN}" >> chrom_recode.log
    echo "  Autosome (1-22) : \${N_AUTO}" >> chrom_recode.log
    echo "  chrX            : \${N_CHRX}" >> chrom_recode.log
    echo "  chrY            : \${N_CHRY}" >> chrom_recode.log
    echo "========================================" >> chrom_recode.log
    echo "SPLIT_CHROM Completed: \$(date)" >> chrom_recode.log
    """
}