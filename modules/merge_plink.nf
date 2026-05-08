process MERGE_PLINK {
    publishDir "${params.outdir}/plink_merged", mode: 'copy'

    input:
    path plink_files   // collected flat list of all per-plate bed/bim/fam files

    output:
    tuple path("cohort.bed"), path("cohort.bim"), path("cohort.fam"), emit: merged

    script:
    """
    # ── Build merge list ───────────────────────────────────────────────────────
    ls *.bed | sed 's/.bed//' | sort > all_prefixes.txt
    head -1 all_prefixes.txt  > base.txt
    tail -n +2 all_prefixes.txt > merge_list.txt

    BASE=\$(cat base.txt)

    if [ -s merge_list.txt ]; then

        # ── Pass 1: attempt merge ──────────────────────────────────────────────
        # Disable bash -e around this call so exit code 3 (3+ alleles) doesn't
        # abort the script before we can inspect and handle it.
        set +e
        plink \\
            --bfile \$BASE \\
            --merge-list merge_list.txt \\
            --make-bed \\
            --allow-extra-chr \\
            --out cohort_tmp
        MERGE1_EXIT=\$?
        set -e

        if [ "\$MERGE1_EXIT" -eq 3 ]; then

            if [ ! -f cohort_tmp-merge.missnp ]; then
                echo "[MERGE_PLINK] ERROR: PLINK exited 3 but cohort_tmp-merge.missnp was not produced." >&2
                exit 1
            fi

            NFLIP=\$(wc -l < cohort_tmp-merge.missnp)
            echo "[MERGE_PLINK] Pass 1: \${NFLIP} variant(s) with 3+ alleles — attempting strand flip." >&2
            echo "[MERGE_PLINK] Affected variant(s):" >&2
            cat cohort_tmp-merge.missnp >&2

            # ── Flip step ─────────────────────────────────────────────────────
            plink \\
                --bfile \$BASE \\
                --flip cohort_tmp-merge.missnp \\
                --make-bed \\
                --allow-extra-chr \\
                --out \${BASE}_flipped

            # ── Pass 2: retry with flipped base ───────────────────────────────
            set +e
            plink \\
                --bfile \${BASE}_flipped \\
                --merge-list merge_list.txt \\
                --make-bed \\
                --allow-extra-chr \\
                --out cohort_tmp
            MERGE2_EXIT=\$?
            set -e

            if [ "\$MERGE2_EXIT" -eq 3 ]; then

                if [ ! -f cohort_tmp-merge.missnp ]; then
                    echo "[MERGE_PLINK] ERROR: Pass 2 exited 3 but cohort_tmp-merge.missnp was not produced." >&2
                    exit 1
                fi

                NREMAIN=\$(wc -l < cohort_tmp-merge.missnp)
                echo "[MERGE_PLINK] WARNING: \${NREMAIN} variant(s) remain multiallelic after strand flip — excluding." >&2
                echo "[MERGE_PLINK] Excluded variant(s):" >&2
                cat cohort_tmp-merge.missnp >&2

                # Save before any subsequent plink call overwrites it
                cp cohort_tmp-merge.missnp multiallelic_exclude.txt

                # Exclude from flipped base
                plink \\
                    --bfile \${BASE}_flipped \\
                    --exclude multiallelic_exclude.txt \\
                    --make-bed \\
                    --allow-extra-chr \\
                    --out \${BASE}_flipped_excl

                # Exclude from every other plate in the merge list
                while read PREFIX; do
                    plink \\
                        --bfile \$PREFIX \\
                        --exclude multiallelic_exclude.txt \\
                        --make-bed \\
                        --allow-extra-chr \\
                        --out \${PREFIX}_excl
                done < merge_list.txt

                # Build new merge list from the cleaned plates
                ls *_excl.bed | sed 's/.bed//' | grep -v "^\${BASE}_flipped_excl\$" | sort > merge_list_excl.txt

                # ── Pass 3: final merge with exclusions ────────────────────────
                set +e
                plink \\
                    --bfile \${BASE}_flipped_excl \\
                    --merge-list merge_list_excl.txt \\
                    --make-bed \\
                    --allow-extra-chr \\
                    --out cohort_tmp
                MERGE3_EXIT=\$?
                set -e

                if [ "\$MERGE3_EXIT" -ne 0 ]; then
                    echo "[MERGE_PLINK] ERROR: Pass 3 (post-exclusion merge) failed with exit code \${MERGE3_EXIT}." >&2
                    exit \$MERGE3_EXIT
                fi

            elif [ "\$MERGE2_EXIT" -ne 0 ]; then
                echo "[MERGE_PLINK] ERROR: Pass 2 merge failed with exit code \${MERGE2_EXIT}." >&2
                exit \$MERGE2_EXIT
            fi

        elif [ "\$MERGE1_EXIT" -ne 0 ]; then
            echo "[MERGE_PLINK] ERROR: Pass 1 merge failed with exit code \${MERGE1_EXIT}." >&2
            exit \$MERGE1_EXIT
        fi

        # ── Final QC: genotype missingness filter ──────────────────────────────
        plink \\
            --bfile cohort_tmp \\
            --geno 0.05 \\
            --make-bed \\
            --allow-extra-chr \\
            --out cohort

    else
        # Only one plate — just rename
        cp \${BASE}.bed cohort.bed
        cp \${BASE}.bim cohort.bim
        cp \${BASE}.fam cohort.fam
    fi
    """
}