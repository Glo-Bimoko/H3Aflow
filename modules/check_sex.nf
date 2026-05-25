process CHECK_SEX {
    tag "cohort"
    publishDir "${params.outdir}/qc/check_sex", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)
    path sex_info

    output:
    path "sexcheck.txt",              emit: sexcheck
    path "sex_discordance.txt",       emit: discordance
    path "sex_plot.png",              emit: plot
    path "sex_check.log",             emit: log
    path "sexcheck_multimind.tsv",    emit: multimind
    path "sexcheck_plate_report.tsv", emit: plate_report

    script:
    def prefix      = bed.baseName
    def bin_dir     = "${projectDir}/bin"
    def sex_info_f  = "${sex_info}"

    def f_lo_male       = params.f_lo_male       ?: 0.8
    def f_hi_female     = params.f_hi_female     ?: 0.2
    def x_maf_male      = params.x_maf_male      ?: 0.01
    def x_maf_female    = params.x_maf_female    ?: 0.01
    def x_miss_male     = params.x_miss_male     ?: 0.8
    def x_miss_female   = params.x_miss_female   ?: 0.05
    def x_diff_miss     = params.x_diff_miss     ?: 0.6
    def x_hwe_female    = params.x_hwe_female    ?: 1e-4
    def x_fisher_p      = params.x_fisher_p      ?: 1e-4
    def mind_thresholds = params.sex_mind_thresholds ?: "0.01,0.03,0.05,0.1"

    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG=sex_check.log
    exec > >(tee -a "\$LOG") 2>&1

    echo "========================================"
    echo "Sex Check Started: \$(date)"
    echo "========================================"
    echo "Parameters:"
    echo "  f_lo_male       : ${f_lo_male}"
    echo "  f_hi_female     : ${f_hi_female}"
    echo "  x_maf_male      : ${x_maf_male}"
    echo "  x_maf_female    : ${x_maf_female}"
    echo "  x_miss_male     : ${x_miss_male}"
    echo "  x_miss_female   : ${x_miss_female}"
    echo "  x_diff_miss     : ${x_diff_miss}"
    echo "  x_hwe_female    : ${x_hwe_female}"
    echo "  x_fisher_p      : ${x_fisher_p}"
    echo "  mind_thresholds : ${mind_thresholds}"
    echo ""

    write_stubs_and_exit() {
        python3 ${bin_dir}/write_sexcheck_stubs.py --reason "\$1"
        python3 ${bin_dir}/write_placeholder_plot.py "\$1"
        echo "Sex Check aborted: \$1"
        exit 0
    }

    N_VAR=\$(wc -l < ${prefix}.bim 2>/dev/null || echo 0)
    echo "chrX variants available (input): \${N_VAR}"
    if [ "\${N_VAR}" -eq 0 ]; then
        write_stubs_and_exit "No chrX variants after QC"
    fi

    echo "--- Step 1: Recode sex_info ---"
    python3 ${bin_dir}/make_sex_update.py ${sex_info_f} plink_sex_update.txt
    echo "  Sex entries to update: \$(wc -l < plink_sex_update.txt)"

    echo "--- Step 2: Rename chr + --split-x b37 ---"
    CHR_NAME=\$(head -1 ${prefix}.bim | cut -f1)
    echo "  Chromosome label in input BIM: '\${CHR_NAME}'"
    if [ "\${CHR_NAME}" = "23" ]; then
        python3 ${bin_dir}/rename_chr23.py ${prefix}.bim chrX_renamed.bim
        ln -sf ${prefix}.bed chrX_renamed.bed
        ln -sf ${prefix}.fam chrX_renamed.fam
        SPLITX_INPUT=chrX_renamed
    else
        SPLITX_INPUT=${prefix}
    fi
    plink --bfile \${SPLITX_INPUT} --split-x b37 no-fail --make-bed --out chrX_splitx
    N_SPLITX=\$(wc -l < chrX_splitx.bim)
    echo "  chrX variants after split-x: \${N_SPLITX} (PAR removed: \$((N_VAR - N_SPLITX)))"

    echo "--- Step 3: Write sex into FAM ---"
    plink --bfile chrX_splitx --update-sex plink_sex_update.txt --make-bed --out chrX_sexed
    N_MALES=\$(awk '\$5==1' chrX_sexed.fam | wc -l)
    N_FEMALES=\$(awk '\$5==2' chrX_sexed.fam | wc -l)
    echo "  Males: \${N_MALES}  Females: \${N_FEMALES}"
    if [ "\${N_MALES}" -eq 0 ] && [ "\${N_FEMALES}" -eq 0 ]; then
        write_stubs_and_exit "Sex update failed — FID mismatch (see sex_check.log)"
    fi

    echo "--- Step 4a: --set-hh-missing ---"
    plink --bfile chrX_sexed --set-hh-missing --make-bed --out chrX_hhmissing

    echo "--- Step 4b: Male freq + miss (post-hh-missing) ---"
    awk '\$5==1 {print \$1, \$2}' chrX_hhmissing.fam > males.txt
    if [ "\$(wc -l < males.txt)" -gt 0 ]; then
        plink --bfile chrX_hhmissing --keep males.txt --freq --missing --allow-no-sex --out chrX_male || true
    fi

    echo "--- Step 4c: Female freq + miss + HWE (pre-hh-missing) ---"
    awk '\$5==2 {print \$1, \$2}' chrX_sexed.fam > females.txt
    if [ "\$(wc -l < females.txt)" -gt 0 ]; then
        plink --bfile chrX_sexed --keep females.txt --freq --missing --hardy --allow-no-sex --out chrX_female || true
    fi

    echo "--- Step 5: chrX SNP QC ---"
    python3 ${bin_dir}/filter_chrx_snps.py \\
        --male_prefix chrX_male \\
        --female_prefix chrX_female \\
        --bim chrX_hhmissing.bim \\
        --out_summary chrX_snp_qc_summary.csv \\
        --out_snps chrX_snps_for_sexcheck.in \\
        --x_maf_male ${x_maf_male} \\
        --x_maf_female ${x_maf_female} \\
        --x_miss_male ${x_miss_male} \\
        --x_miss_female ${x_miss_female} \\
        --x_diff_miss ${x_diff_miss} \\
        --x_hwe_female ${x_hwe_female} \\
        --x_fisher_p ${x_fisher_p}

    echo "--- Step 6: Extract QC-passing variants ---"
    N_PASS=\$(wc -l < chrX_snps_for_sexcheck.in)
    if [ "\${N_PASS}" -eq 0 ]; then
        echo "  WARNING: 0 SNPs passed — falling back to full chrX_hhmissing"
        FINAL_BED=chrX_hhmissing
    else
        plink --bfile chrX_hhmissing --extract chrX_snps_for_sexcheck.in --make-bed --out chrX_for_sexcheck
        FINAL_BED=chrX_for_sexcheck
    fi
    echo "  Final chrX variant count: \$(wc -l < "\${FINAL_BED}.bim")"

    echo "--- Step 7: PLINK --check-sex (base run) ---"
    plink --bfile "\${FINAL_BED}" --check-sex ${f_hi_female} ${f_lo_male} --out sexcheck || true
    if [ ! -f sexcheck.sexcheck ]; then
        printf 'FID IID PEDSEX SNPSEX STATUS F\\n' > sexcheck.sexcheck
    fi
    echo "  PLINK STATUS=OK: \$(awk 'NR>1 && \$5==\"OK\"' sexcheck.sexcheck | wc -l)"
    echo "  PLINK STATUS=PROBLEM: \$(awk 'NR>1 && \$5==\"PROBLEM\"' sexcheck.sexcheck | wc -l)"

    echo "--- Step 7b: Multi-missingness --check-sex ---"
    python3 ${bin_dir}/sexcheck_multimind.py \\
        --bfile "\${FINAL_BED}" \\
        --out sexcheck_multimind.tsv \\
        --f_hi_female ${f_hi_female} \\
        --f_lo_male ${f_lo_male} \\
        --mind_thresholds ${mind_thresholds}

    echo "--- Step 7c: Per-plate discordance report ---"
    python3 ${bin_dir}/sexcheck_plate_report.py \\
        --sexcheck sexcheck.sexcheck \\
        --multimind sexcheck_multimind.tsv \\
        --sex_info ${sex_info_f} \\
        --out sexcheck_plate_report.tsv

    echo "--- Step 8: Annotate with collected sex ---"
    python3 ${bin_dir}/annotate_sex_check.py \\
        --sexcheck sexcheck.sexcheck \\
        --sex_info ${sex_info_f} \\
        --f_lo_male ${f_lo_male} \\
        --f_hi_female ${f_hi_female} \\
        --out_annot sexcheck.txt \\
        --out_discord sex_discordance.txt || python3 ${bin_dir}/write_sexcheck_stubs.py --reason "annotate_sex_check.py failed"

    echo "--- Step 9: Generate sex check plot ---"
    python3 ${bin_dir}/compare_sex.py \\
        --plink_sex sexcheck.sexcheck \\
        --collected_sex ${sex_info_f} \\
        --f_lo_male ${f_lo_male} \\
        --f_hi_female ${f_hi_female} \\
        --out sex_discordance_plot.tsv \\
        --plot sex_plot.png || python3 ${bin_dir}/write_placeholder_plot.py "Sex check plot unavailable (see sex_check.log)"

    echo "========================================"
    echo "Sex Check Summary"
    echo "  Collected≠inferred discordant : \$(awk 'NR>1' sex_discordance.txt | wc -l)"
    echo "  Hard discordant               : \$(awk -F'\\t' 'NR>1 && \$NF==\"HARD_DISCORDANT\"' sexcheck_multimind.tsv | wc -l)"
    echo "  Variable discordant           : \$(awk -F'\\t' 'NR>1 && \$NF==\"DISCORDANT_VARIABLE\"' sexcheck_multimind.tsv | wc -l)"
    echo "Sex Check Completed: \$(date)"
    echo "========================================"
    """
}
