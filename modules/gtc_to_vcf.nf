process GTC_TO_VCF {
    tag "$plate"
    publishDir "${params.outdir}/bcf", mode: 'copy', pattern: "*.bcf*"
    publishDir "${params.outdir}/tsv", mode: 'copy', pattern: "*.tsv"
   
    input:
    tuple val(plate), path(gtc_files)
    val bpm
    val egt
    val fasta

    output:
    tuple val(plate), path("${plate}.bcf"),     emit: bcf
    tuple val(plate), path("${plate}.bcf.csi"), emit: bcf_index
    tuple val(plate), path("${plate}.tsv"),     emit: tsv

    script:
    // Sanitise plate name for use as filename (replace spaces/special chars)
    def safe_plate = plate.replaceAll(/[^A-Za-z0-9_\-]/, '_')
    """
    # Ensure the BCFTOOLS_PLUGINS environment variable is set for the plugin to be found
    export BCFTOOLS_PLUGINS=${projectDir}/bin/plugins/gtc2vcf

    ls -1 *.gtc > gtc_list.txt

    python ${projectDir}/bin/convert_gtc2vcf.py \
        --bpm       "${bpm}"       \
        --egt       "${egt}"       \
        --gtcs      gtc_list.txt   \
        --fasta     "${fasta}"     \
        --outprefix "${safe_plate}"

    # Rename outputs back to the original plate name if different
    if [ "${safe_plate}" != "${plate}" ]; then
        mv "${safe_plate}.bcf" "${plate}.bcf"     2>/dev/null || true
        mv "${safe_plate}.bcf.csi" "${plate}.bcf.csi" 2>/dev/null || true
        mv "${safe_plate}.tsv" "${plate}.tsv"     2>/dev/null || true
    fi
    """
}
