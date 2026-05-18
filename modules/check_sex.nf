process CHECK_SEX {
    tag "cohort"
    publishDir "${params.outdir}/qc/check_sex", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)
    path sex_info

    output:
    path "sexcheck.txt",        emit: sexcheck
    path "sex_discordance.txt", emit: discordance
    path "sex_plot.png",        emit: plot
    path "sex_check.log",       emit: log

    script:
    def prefix         = bed.baseName
    def annot_script   = "${projectDir}/bin/annotate_sex_check.py"
    def compare_script = "${projectDir}/bin/compare_sex.py"
    def sex_update_py  = "${projectDir}/bin/make_sex_update.py"
    def placeholder_py = "${projectDir}/bin/write_placeholder_plot.py"
    def sex_info_path  = "${sex_info}"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    echo "========================================" > sex_check.log
    echo "Sex Check Started: \$(date)" >> sex_check.log
    echo "========================================" >> sex_check.log

    # ── Preflight: count chrX variants ───────────────────────────────────────
    N_VAR=\$(wc -l < ${prefix}.bim 2>/dev/null || echo 0)
    echo "chrX variants available: \${N_VAR}" >> sex_check.log

    if [ "\${N_VAR}" -eq 0 ]; then
        echo "ERROR: No chrX variants found. Sex check cannot run." >> sex_check.log
        printf 'FID\tIID\tPEDSEX\tSNPSEX\tSTATUS\tF\tCOLLECTED_SEX\tINFERRED_SEX\tDISCORDANT\n' \
            > sexcheck.txt
        cp sexcheck.txt sex_discordance.txt
        python3 ${placeholder_py} "Sex check could not run (no chrX variants after QC)"
        exit 0
    fi

    # ── Step 1: Recode sex_info to PLINK format ───────────────────────────────
    # make_sex_update.py now writes FID=IID (matching --double-id convention).
    echo "Step 1: Recoding sex_info to PLINK format..." >> sex_check.log
    python3 ${sex_update_py} ${sex_info_path} plink_sex_update.txt
    N_SEX=\$(wc -l < plink_sex_update.txt)
    echo "  Sex entries to update: \${N_SEX}" >> sex_check.log

    # Diagnostic: show first few lines of update file and FAM so FID mismatch is obvious in logs
    echo "  plink_sex_update.txt (first 3 lines):" >> sex_check.log
    head -3 plink_sex_update.txt >> sex_check.log 2>&1 || true
    echo "  ${prefix}.fam (first 3 lines):" >> sex_check.log
    head -3 ${prefix}.fam >> sex_check.log 2>&1 || true

    # ── Step 2: Split PAR from chrX ───────────────────────────────────────────
    # PAR1/PAR2 are diploid in both sexes; they bias male F toward 0 if kept.
    # --split-x b37 moves them to a separate XY contig.
    # Note: VCF_TO_PLINK already runs --split-x, but we repeat it here as a
    # safety net in case the merged cohort dataset lost the XY contig labeling.
    echo "Step 2: Splitting PAR from chrX (--split-x b37)..." >> sex_check.log
    CHR_NAME=\$(head -1 ${prefix}.bim | cut -f1)
    echo "  Chromosome name in BIM: \${CHR_NAME}" >> sex_check.log

    if [ "\${CHR_NAME}" = "23" ]; then
        echo "  Numeric chr naming — remapping 23->X for split-x..." >> sex_check.log
        python3 ${projectDir}/bin/rename_chr23.py ${prefix}.bim chrX_renamed.bim
        ln -sf ${prefix}.bed chrX_renamed.bed
        ln -sf ${prefix}.fam chrX_renamed.fam
        plink \\
            --bfile chrX_renamed \\
            --split-x b37 no-fail \\
            --make-bed \\
            --out chrX_splitx >> sex_check.log 2>&1
    else
        plink \\
            --bfile ${prefix} \\
            --split-x b37 no-fail \\
            --make-bed \\
            --out chrX_splitx >> sex_check.log 2>&1
    fi

    N_VAR_SPLITX=\$(wc -l < chrX_splitx.bim 2>/dev/null || echo 0)
    echo "  chrX variants after split-x: \${N_VAR_SPLITX}" >> sex_check.log

    # ── Step 3: Write sex into FAM ────────────────────────────────────────────
    # --update-sex must be a separate PLINK call from --check-sex in PLINK 1.9.
    echo "Step 3: Writing sex into FAM..." >> sex_check.log
    plink \\
        --bfile chrX_splitx \\
        --update-sex plink_sex_update.txt \\
        --make-bed \\
        --out chrX_sexed >> sex_check.log 2>&1

    # FAM sex is column 5 (not 3 — col 3 is paternal ID)
    N_MALES=\$(awk '\$5==1' chrX_sexed.fam | wc -l)
    N_FEMALES=\$(awk '\$5==2' chrX_sexed.fam | wc -l)
    echo "  Males in FAM   : \${N_MALES}" >> sex_check.log
    echo "  Females in FAM : \${N_FEMALES}" >> sex_check.log

    # Guard: if sex update didn't take, the rest of the sex check is meaningless
    if [ "\${N_MALES}" -eq 0 ] && [ "\${N_FEMALES}" -eq 0 ]; then
        echo "" >> sex_check.log
        echo "ERROR: Sex update failed — 0 males and 0 females in FAM after --update-sex." >> sex_check.log
        echo "  This means FID values in plink_sex_update.txt do not match chrX_splitx.fam." >> sex_check.log
        echo "  Check that VCF_TO_PLINK uses --double-id and make_sex_update.py writes FID=IID." >> sex_check.log
        echo "  First 5 lines of chrX_splitx.fam:" >> sex_check.log
        head -5 chrX_splitx.fam >> sex_check.log 2>&1 || true
        echo "  First 5 lines of plink_sex_update.txt:" >> sex_check.log
        head -5 plink_sex_update.txt >> sex_check.log 2>&1 || true
        python3 ${placeholder_py} "Sex update failed — FID mismatch (see sex_check.log)"
        # Write minimal valid output files so the pipeline doesn't crash downstream
        printf 'FID\tIID\tPEDSEX\tSNPSEX\tSTATUS\tF\tCOLLECTED_SEX\tINFERRED_SEX\tDISCORDANT\n' > sexcheck.txt
        printf 'FID\tIID\tCOLLECTED_SEX\tINFERRED_SEX\tF\tSTATUS\n' > sex_discordance.txt
        exit 0
    fi

    # ── Step 4: Identify het-haploid variants (preflight) ────────────────────
    echo "Step 4: Identifying het-haploid variants (preflight)..." >> sex_check.log
    plink \\
        --bfile chrX_sexed \\
        --check-sex \\
        --out sexcheck_preflight >> sex_check.log 2>&1 || true

    N_HH=0
    HH_FILE=""
    for candidate in chrX_sexed.hh sexcheck_preflight.hh; do
        if [ -f "\${candidate}" ]; then
            HH_FILE="\${candidate}"
            N_HH=\$(wc -l < "\${HH_FILE}")
            break
        fi
    done
    echo "  Het-haploid genotype lines found: \${N_HH}" >> sex_check.log

    FINAL_BED=chrX_sexed

    if [ "\${N_HH}" -gt 0 ]; then
        awk '{print \$3}' "\${HH_FILE}" | sort -u > hh_snps.txt
        N_HH_SNPS=\$(wc -l < hh_snps.txt)
        echo "  Unique het-haploid variant IDs: \${N_HH_SNPS}" >> sex_check.log

        plink \\
            --bfile chrX_sexed \\
            --exclude hh_snps.txt \\
            --make-bed \\
            --out chrX_sexed_nohh >> sex_check.log 2>&1

        N_VAR_NOHH=\$(wc -l < chrX_sexed_nohh.bim 2>/dev/null || echo 0)
        echo "  chrX variants after het-haploid exclusion: \${N_VAR_NOHH}" >> sex_check.log

        if [ "\${N_VAR_NOHH}" -eq 0 ]; then
            echo "" >> sex_check.log
            echo "WARNING: Zero chrX variants remain after het-haploid exclusion." >> sex_check.log
            echo "  All chrX variants were het-haploid — haploid recoding never ran." >> sex_check.log
            echo "  Root cause: convert_gtc2vcf.py Step 3 (+setGT) did not execute." >> sex_check.log
            echo "  Check that --sex-info was passed to GTC_TO_VCF." >> sex_check.log
            echo "  Using preflight sexcheck as best-available fallback." >> sex_check.log
            if [ -f sexcheck_preflight.sexcheck ]; then
                cp sexcheck_preflight.sexcheck sexcheck.sexcheck
            else
                printf 'FID IID PEDSEX SNPSEX STATUS F\n' > sexcheck.sexcheck
            fi
            FINAL_BED=""
        else
            FINAL_BED=chrX_sexed_nohh
        fi
    fi

    # ── Step 5: Run PLINK --check-sex ─────────────────────────────────────────
    if [ -n "\${FINAL_BED}" ]; then
        echo "Step 5: Checking for polymorphic variants before --check-sex..." >> sex_check.log
        plink --bfile "\${FINAL_BED}" --freq --out freq_check >> sex_check.log 2>&1 || true
        N_POLY=0
        if [ -f freq_check.frq ]; then
            N_POLY=\$(awk 'NR>1 && \$5+0>0 && \$5+0<1' freq_check.frq | wc -l)
        fi
        echo "  Polymorphic chrX variants: \${N_POLY}" >> sex_check.log

        if [ "\${N_POLY}" -eq 0 ]; then
            echo "" >> sex_check.log
            echo "WARNING: No polymorphic chrX variants — cannot compute F-statistic." >> sex_check.log
            echo "  Most likely cause: haploid recoding did not run in convert_gtc2vcf.py." >> sex_check.log
            echo "  Using preflight sexcheck as fallback (results unreliable)." >> sex_check.log
            if [ -f sexcheck_preflight.sexcheck ]; then
                cp sexcheck_preflight.sexcheck sexcheck.sexcheck
            else
                printf 'FID IID PEDSEX SNPSEX STATUS F\n' > sexcheck.sexcheck
            fi
        else
            echo "Step 5: Running PLINK --check-sex on \${FINAL_BED}..." >> sex_check.log
            plink \\
                --bfile "\${FINAL_BED}" \\
                --check-sex \\
                --out sexcheck >> sex_check.log 2>&1

            if [ ! -f sexcheck.sexcheck ]; then
                echo "ERROR: PLINK did not produce sexcheck.sexcheck." >> sex_check.log
                printf 'FID IID PEDSEX SNPSEX STATUS F\n' > sexcheck.sexcheck
            fi
        fi
    fi

    # ── Count results ─────────────────────────────────────────────────────────
    N_PROBLEM=\$(awk 'NR>1 && \$5=="PROBLEM"' sexcheck.sexcheck | wc -l || echo 0)
    N_OK=\$(awk 'NR>1 && \$5=="OK"' sexcheck.sexcheck | wc -l || echo 0)
    N_TOTAL=\$(awk 'NR>1' sexcheck.sexcheck | wc -l || echo 0)
    N_F_NEG1=\$(awk 'NR>1 && \$6=="-1"' sexcheck.sexcheck | wc -l || echo 0)
    echo "  PLINK STATUS=OK     : \${N_OK}" >> sex_check.log
    echo "  PLINK STATUS=PROBLEM: \${N_PROBLEM}" >> sex_check.log

    if [ "\${N_TOTAL}" -gt 0 ] && [ "\${N_F_NEG1}" -eq "\${N_TOTAL}" ]; then
        echo "WARNING: All samples have F=-1 — male chrX still diploid." >> sex_check.log
        echo "  Re-run GTC_TO_VCF with --sex-info passed to convert_gtc2vcf.py." >> sex_check.log
    fi

    # ── Step 6: Annotate with collected sex ───────────────────────────────────
    echo "Step 6: Annotating with collected sex..." >> sex_check.log
    python3 ${annot_script} \\
        --sexcheck    sexcheck.sexcheck \\
        --sex_info    ${sex_info_path} \\
        --out_annot   sexcheck.txt \\
        --out_discord sex_discordance.txt >> sex_check.log 2>&1 || {
        echo "WARNING: annotate_sex_check.py failed — writing empty fallback files." >> sex_check.log
        printf 'FID\tIID\tPEDSEX\tSNPSEX\tSTATUS\tF\tCOLLECTED_SEX\tINFERRED_SEX\tDISCORDANT\n' \
            > sexcheck.txt
        printf 'FID\tIID\tCOLLECTED_SEX\tINFERRED_SEX\tF\tSTATUS\n' > sex_discordance.txt
    }

    # ── Step 7: Generate plot ─────────────────────────────────────────────────
    echo "Step 7: Generating sex check plot..." >> sex_check.log
    python3 ${compare_script} \\
        --plink_sex     sexcheck.sexcheck \\
        --collected_sex ${sex_info_path} \\
        --out           sex_discordance_plot.tsv \\
        --plot          sex_plot.png >> sex_check.log 2>&1 || {
        echo "WARNING: compare_sex.py failed — generating placeholder plot." >> sex_check.log
        python3 ${placeholder_py} "Sex check plot unavailable (see sex_check.log)"
    }

    # ── Summary ───────────────────────────────────────────────────────────────
    N_DISCORD=\$(awk 'NR>1' sex_discordance.txt 2>/dev/null | wc -l || echo 0)
    N_VALID_F=\$(awk 'NR>1 && \$6!="nan" && \$6!="NA" && \$6!="" && \$6!="-1"' \\
                    sexcheck.sexcheck 2>/dev/null | wc -l || echo 0)

    echo "" >> sex_check.log
    echo "========================================" >> sex_check.log
    echo "Sex Check Summary:" >> sex_check.log
    echo "  chrX variants (raw)          : \${N_VAR}" >> sex_check.log
    echo "  chrX variants (after split-x): \${N_VAR_SPLITX:-n/a}" >> sex_check.log
    echo "  Het-haploid genotype lines   : \${N_HH}" >> sex_check.log
    echo "  Samples sex updated          : \${N_SEX}" >> sex_check.log
    echo "  Males in FAM                 : \${N_MALES}" >> sex_check.log
    echo "  Females in FAM               : \${N_FEMALES}" >> sex_check.log
    echo "  PLINK OK                     : \${N_OK}" >> sex_check.log
    echo "  PLINK PROBLEM                : \${N_PROBLEM}" >> sex_check.log
    echo "  Samples with valid F-stat    : \${N_VALID_F}" >> sex_check.log
    echo "  True discordant              : \${N_DISCORD}  (collected vs inferred)" >> sex_check.log
    echo "========================================" >> sex_check.log
    echo "Sex Check Completed: \$(date)" >> sex_check.log

    if [ "\${N_VALID_F}" -eq 0 ] && [ "\${N_TOTAL}" -gt 0 ]; then
        echo "" >> sex_check.log
        echo "ACTION REQUIRED: No valid F-statistics produced." >> sex_check.log
        echo "  Sex check results are unreliable for this run." >> sex_check.log
        echo "  Fix: re-run from GTC_TO_VCF after deploying the updated" >> sex_check.log
        echo "  convert_gtc2vcf.py, which recodes male chrX genotypes to haploid." >> sex_check.log
    fi
    """
}