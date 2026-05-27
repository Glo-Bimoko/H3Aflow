process SAMPLE_QC {
    publishDir "${params.outdir}/qc/sample_qc", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    path "sample_qc.imiss",        emit: imiss
    path "sample_qc.het",          emit: het
    path "sample_qc_stats.tsv",    emit: qc_stats
    path "samples_pass_qc.txt",    emit: keep_list
    path "samples_fail_qc.txt",    emit: fail_list

    script:
    """
    PREFIX=\$(basename ${bed} .bed)

    # ── 1. Missingness per sample ──────────────────────────────────────────
    plink \
        --bfile \${PREFIX} \
        --missing \
        --allow-extra-chr \
        --out sample_qc

    # ── 2. Heterozygosity ─────────────────────────────────────────────────
    plink \
        --bfile \${PREFIX} \
        --het \
        --allow-extra-chr \
        --out sample_qc

    # ── 3. Compute combined stats + flag failures ──────────────────────────
    python ${projectDir}/bin/compute_sample_qc.py \
        --imiss  sample_qc.imiss \
        --het    sample_qc.het   \
        --mind   ${params.mind}  \
        --het_sd ${params.het_sd}\
        --out_stats  sample_qc_stats.tsv \
        --out_pass   samples_pass_qc.txt \
        --out_fail   samples_fail_qc.txt
    """
}
