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
    // Sanitise plate name for use as filename (replace spaces/special chars)
    def safe_plate = plate.replaceAll(/[^A-Za-z0-9_\-]/, '_')
    """
    # ── Plugin paths ──────────────────────────────────────────────────────────
    # APPEND the gtc2vcf plugin dir to whatever BCFTOOLS_PLUGINS is already set
    # to by the Nextflow profile (e.g. the standard bcftools libexec dir added
    # via env.BCFTOOLS_PLUGINS in nextflow.config).  Using += avoids clobbering
    # the standard plugin dir that +setGT and +fill-tags live in.
    export BCFTOOLS_PLUGINS=\${BCFTOOLS_PLUGINS:+\${BCFTOOLS_PLUGINS}:}${projectDir}/bin/plugins/gtc2vcf

    # Forward the explicit libexec override so Python's build_plugin_env() can
    # pick it up regardless of how the Nextflow profile propagates env vars.
    export BCFTOOLS_LIBEXEC=\${BCFTOOLS_LIBEXEC:-}

    ls -1 *.gtc > gtc_list.txt

    python ${projectDir}/bin/convert_gtc2vcf.py \
        --bpm       "${bpm}"        \
        --egt       "${egt}"        \
        --gtcs      gtc_list.txt    \
        --fasta-ref "${fasta}"      \
        --outprefix "${safe_plate}" \
        --sex-info  "${sex_info}"

    # Rename outputs back to the original plate name if different
    if [ "${safe_plate}" != "${plate}" ]; then
        mv "${safe_plate}.bcf"     "${plate}.bcf"     2>/dev/null || true
        mv "${safe_plate}.bcf.csi" "${plate}.bcf.csi" 2>/dev/null || true
        mv "${safe_plate}.tsv"     "${plate}.tsv"     2>/dev/null || true
    fi
    """
}