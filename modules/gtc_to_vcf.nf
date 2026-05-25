process GTC_TO_VCF {
    tag "$plate"
    publishDir "${params.outdir}/bcf", mode: 'copy', pattern: "*.bcf*"
    publishDir "${params.outdir}/tsv", mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(plate), path(gtc_files)
    val bpm
    val egt
    val fasta
    path sex_info

    output:
    tuple val(plate), path("${plate}.bcf"),     emit: bcf
    tuple val(plate), path("${plate}.bcf.csi"), emit: bcf_index
    tuple val(plate), path("${plate}.tsv"),     emit: tsv

    script:
    def safe_plate = plate.replaceAll(/[^A-Za-z0-9_\-]/, '_')
    """
    set -euo pipefail

    # ── 0. Plugin path ────────────────────────────────────────────────────────
    # Only +gtc2vcf is needed here — no +setGT.
    # We do NOT recode male chrX to haploid at the BCF stage.
    # Rationale: bcftools +setGT has proven fragile across versions.
    # PLINK handles male chrX ploidy correctly when sex is written into the FAM
    # via --update-sex (done in CHECK_SEX). The het-haploid variants that PLINK
    # flags are excluded by the preflight step in check_sex.nf before --check-sex
    # runs, so the sex check result is valid without BCF-level recoding.
    export BCFTOOLS_PLUGINS="${projectDir}/bin/plugins/gtc2vcf\${BCFTOOLS_PLUGINS:+:\${BCFTOOLS_PLUGINS}}"
    echo "[GTC_TO_VCF] BCFTOOLS_PLUGINS=\${BCFTOOLS_PLUGINS}"

    # ── 1. GTC → sorted, normalised BCF ──────────────────────────────────────
    # Following Ayoub Ksouri's proven approach:
    #   +gtc2vcf | sort | norm → indexed BCF
    # --extra writes XY intensities to TSV during conversion (no separate step).
    # bcftools norm -c x: if REF mismatches the FASTA, set GT to missing rather
    # than aborting — keeps the pipeline robust to manifest/reference mismatches.
    echo "[GTC_TO_VCF] Step 1: GTC -> sorted, normalised BCF"
    ls -1 *.gtc > gtc_list.txt
    echo "[GTC_TO_VCF] GTCs: \$(wc -l < gtc_list.txt)"

    bcftools +gtc2vcf \\
        --no-version \\
        --output-type u \\
        --bpm       "${bpm}" \\
        --egt       "${egt}" \\
        --fasta-ref "${fasta}" \\
        --gtcs      gtc_list.txt \\
        --extra     "${safe_plate}.tsv" | \\
    bcftools sort \\
        --output-type u \\
        --temp-dir  . | \\
    bcftools norm \\
        --no-version \\
        --output-type b \\
        --output    "${safe_plate}.bcf" \\
        --fasta-ref "${fasta}" \\
        --check-ref x \\
        --write-index

    echo "[GTC_TO_VCF] BCF written: ${safe_plate}.bcf"

    # ── 2. Rename if plate name had special characters ────────────────────────
    if [ "${safe_plate}" != "${plate}" ]; then
        mv "${safe_plate}.bcf"     "${plate}.bcf"     2>/dev/null || true
        mv "${safe_plate}.bcf.csi" "${plate}.bcf.csi" 2>/dev/null || true
        mv "${safe_plate}.tsv"     "${plate}.tsv"     2>/dev/null || true
    fi

    # ── 3. Verify BCF and TSV ─────────────────────────────────────────────────
    N_SAMPLES=\$(bcftools query --list-samples "${plate}.bcf" | wc -l)
    N_VARS=\$(bcftools view -H "${plate}.bcf" | wc -l)
    echo "[GTC_TO_VCF] Samples: \${N_SAMPLES}  Variants: \${N_VARS}"

    if [ ! -f "${plate}.tsv" ] || [ ! -s "${plate}.tsv" ]; then
        echo "[GTC_TO_VCF] WARNING: TSV missing or empty — XY intensity QC will be skipped" >&2
        echo -e "SAMPLE_ID\\tCHR\\tPOS\\tREF\\tALT\\tNORMX\\tNORMY" > "${plate}.tsv"
    fi

    rm -f gtc_list.txt
    echo "[GTC_TO_VCF] Done — plate ${plate}"
    """
}