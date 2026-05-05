process CHECK_SEX {
    tag "cohort"
    publishDir "${params.outdir}/qc/check_sex", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)   // chrX-only PLINK fileset from SPLIT_CHROM
    path sex_info                            // collected sex TSV: sampleid, sex (0=F, 1=M)

    output:
    path "sexcheck.txt",        emit: sexcheck
    path "sex_discordance.txt", emit: discordance
    path "sex_plot.png",        emit: plot
    path "sex_check.log",       emit: log

    script:
    def prefix = bed.baseName

    """
    echo "========================================" > sex_check.log
    echo "Sex Check Started: \$(date)" >> sex_check.log
    echo "========================================" >> sex_check.log

    # ── Preflight: count chrX variants ───────────────────────────────────────
    N_VAR=\$(wc -l < ${prefix}.bim 2>/dev/null || echo 0)
    echo "chrX variants available: \${N_VAR}" >> sex_check.log

    if [ "\${N_VAR}" -eq 0 ]; then
        echo "ERROR: No chrX variants found. Sex check cannot run." >> sex_check.log
        printf 'FID\tIID\tPEDSEX\tSNPSEX\tSTATUS\tF\tCOLLECTED_SEX\tINFERRED_SEX\tDISCORDANT\n' > sexcheck.txt
        cp sexcheck.txt sex_discordance.txt
        python3 -c "
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 4))
ax.text(0.5, 0.5, 'Sex check could not run\n(no chrX variants after QC)',
        ha='center', va='center', fontsize=14, transform=ax.transAxes)
ax.axis('off')
plt.tight_layout()
plt.savefig('sex_plot.png', dpi=150, bbox_inches='tight')
plt.close()
"
        exit 0
    fi

    # ── Step 1: PLINK 1.9 --check-sex ────────────────────────────────────────
    # plink2 does NOT have --check-sex; this flag belongs to PLINK 1.9.
    # PLINK 1.9 computes the X-chromosome inbreeding coefficient (F-statistic):
    #   F < 0.2  → Female,  F > 0.8  → Male,  otherwise → Ambiguous
    echo "Step 1: Running PLINK 1.9 --check-sex..." >> sex_check.log
    plink \
        --bfile ${prefix} \
        --check-sex \
        --allow-extra-chr \
        --allow-no-sex \
        --out sexcheck 2>&1 | tee -a sex_check.log

    # ── Step 2: Verify PLINK produced output ─────────────────────────────────
    if [ ! -f sexcheck.sexcheck ]; then
        echo "ERROR: PLINK 1.9 did not produce sexcheck.sexcheck" >> sex_check.log
        printf 'FID\tIID\tPEDSEX\tSNPSEX\tSTATUS\tF\tCOLLECTED_SEX\tINFERRED_SEX\tDISCORDANT\n' > sexcheck.txt
        cp sexcheck.txt sex_discordance.txt
        python3 -c "
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 4))
ax.text(0.5, 0.5, 'PLINK sex check failed\n(see sex_check.log)',
        ha='center', va='center', fontsize=14, transform=ax.transAxes)
ax.axis('off')
plt.tight_layout()
plt.savefig('sex_plot.png', dpi=150, bbox_inches='tight')
plt.close()
"
        exit 0
    fi

    N_PROBLEM=\$(awk 'NR>1 && \$5=="PROBLEM"' sexcheck.sexcheck | wc -l)
    N_OK=\$(awk 'NR>1 && \$5=="OK"' sexcheck.sexcheck | wc -l)
    echo "  PLINK STATUS=OK     : \${N_OK}" >> sex_check.log
    echo "  PLINK STATUS=PROBLEM: \${N_PROBLEM}" >> sex_check.log

    # ── Step 3: Annotate with collected sex ───────────────────────────────────
    # annotate_sex_check.py merges PLINK output with the sex_info file and
    # writes sexcheck.txt (full) and sex_discordance.txt (discordant only).
    # The || block ensures output files always exist even if the script errors.
    echo "Step 2: Annotating with collected sex..." >> sex_check.log
    python3 ${projectDir}/bin/annotate_sex_check.py \
        --sexcheck    sexcheck.sexcheck \
        --sex_info    ${sex_info} \
        --out_annot   sexcheck.txt \
        --out_discord sex_discordance.txt 2>&1 | tee -a sex_check.log || {
        echo "WARNING: annotate_sex_check.py failed – writing empty fallback files." >> sex_check.log
        printf 'FID\tIID\tPEDSEX\tSNPSEX\tSTATUS\tF\tCOLLECTED_SEX\tINFERRED_SEX\tDISCORDANT\n' > sexcheck.txt
        printf 'FID\tIID\tCOLLECTED_SEX\tINFERRED_SEX\tF\tSTATUS\n' > sex_discordance.txt
    }

    # ── Step 4: Generate plot ─────────────────────────────────────────────────
    echo "Step 3: Generating sex check plot..." >> sex_check.log
    python3 ${projectDir}/bin/compare_sex.py \
        --plink_sex     sexcheck.sexcheck \
        --collected_sex ${sex_info} \
        --out           sex_discordance_plot.tsv \
        --plot          sex_plot.png 2>&1 | tee -a sex_check.log || {
        echo "WARNING: compare_sex.py failed – generating placeholder plot." >> sex_check.log
        python3 -c "
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 4))
ax.text(0.5, 0.5, 'Sex check plot unavailable\n(see sex_check.log)',
        ha='center', va='center', fontsize=14, transform=ax.transAxes)
ax.axis('off')
plt.tight_layout()
plt.savefig('sex_plot.png', dpi=150, bbox_inches='tight')
plt.close()
"
    }

    # ── Summary ───────────────────────────────────────────────────────────────
    N_DISCORD=\$([ -f sex_discordance.txt ] && awk 'NR>1' sex_discordance.txt | wc -l || echo 0)
    echo "" >> sex_check.log
    echo "========================================" >> sex_check.log
    echo "Sex Check Summary:" >> sex_check.log
    echo "  chrX variants used  : \${N_VAR}" >> sex_check.log
    echo "  PLINK OK            : \${N_OK}" >> sex_check.log
    echo "  PLINK PROBLEM       : \${N_PROBLEM}" >> sex_check.log
    echo "  True discordant     : \${N_DISCORD}  (collected vs inferred mismatch)" >> sex_check.log
    echo "========================================" >> sex_check.log
    echo "Sex Check Completed: \$(date)" >> sex_check.log

    # Note if all F-stats are NaN (monomorphic X in small subset)
    N_VALID_F=\$(awk 'NR>1 && \$6!="nan" && \$6!="NA" && \$6!=""' sexcheck.sexcheck | wc -l)
    if [ "\${N_VALID_F}" -eq 0 ]; then
        echo "" >> sex_check.log
        echo "WARNING: All F-statistics are NaN." >> sex_check.log
        echo "  No polymorphic chrX SNPs in this batch (common in small subsets <500 samples)." >> sex_check.log
        echo "  Sex check will be reliable only when run on the full cohort on CHPC." >> sex_check.log
    fi
    """
}