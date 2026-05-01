process VCF_TO_PLINK {
    tag "$plate"
    publishDir "${params.outdir}/plink_per_plate", mode: 'copy'

    input:
    tuple val(plate), path(bcf), path(bcf_index)

    output:
    tuple val(plate), path("${plate}.bed"),
                      path("${plate}.bim"),
                      path("${plate}.fam"), emit: plink

    script:
    // Sanitise plate name — spaces and special chars break plink --out
    def safe = plate.replaceAll(/[^A-Za-z0-9_\-]/, '_')
    """
    # Convert BCF → PLINK binary format
    plink \
        --bcf         "${bcf}"   \
        --keep-allele-order       \
        --vcf-idspace-to _        \
        --const-fid 0             \
        --allow-extra-chr         \
        --split-x b37 no-fail     \
        --make-bed                \
        --out        "${safe}"

    # Rename to original plate name if sanitised name differs
    if [ "${safe}" != "${plate}" ]; then
        for ext in bed bim fam log; do
            [ -f "${safe}.\${ext}" ] && mv "${safe}.\${ext}" "${plate}.\${ext}" || true
        done
    fi
    """
}