process MERGE_PLINK {
    publishDir "${params.outdir}/plink_merged", mode: 'copy'

    input:
    path plink_files   // collected list of all per-plate bed/bim/fam files

    output:
    tuple path("cohort.bed"), path("cohort.bim"), path("cohort.fam"), emit: merged

    script:
    """
    # Build a merge-list: one bed prefix per line (excluding the first plate,
    # which is the base dataset passed to --bfile)
    ls *.bed | sed 's/.bed//' | sort > all_prefixes.txt
    head -1 all_prefixes.txt  > base.txt
    tail -n +2 all_prefixes.txt > merge_list.txt

    BASE=\$(cat base.txt)

    if [ -s merge_list.txt ]; then
        plink \
            --bfile \$BASE \
            --merge-list merge_list.txt \
            --make-bed \
            --allow-extra-chr \
            --out cohort_tmp

        # Remove ambiguous (A/T, G/C) SNPs and apply basic variant QC
        plink \
            --bfile cohort_tmp \
            --geno 0.05 \
            --make-bed \
            --allow-extra-chr \
            --out cohort
    else
        # Only one plate – just rename
        cp \${BASE}.bed cohort.bed
        cp \${BASE}.bim cohort.bim
        cp \${BASE}.fam cohort.fam
    fi
    """
}
