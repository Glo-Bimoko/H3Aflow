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
    #
    # Key flags (aligned with collaborator's idat2plink.sh):
    #   --double-id         : sets FID = IID (avoids FID=0 mismatch in --update-sex)
    #   --keep-allele-order : preserves ref/alt orientation from the BCF
    #   --allow-extra-chr   : passes non-standard contigs through without error
    #   --split-x b37       : moves PAR1/PAR2 to XY contig so male X is truly haploid
    #
    # NOTE: --const-fid 0 removed; --double-id sets FID = IID for downstream PLINK.
    plink \
        --bcf              "${bcf}"   \
        --keep-allele-order           \
        --vcf-idspace-to _            \
        --double-id                   \
        --allow-extra-chr             \
        --split-x b37 no-fail         \
        --make-bed                    \
        --out              "${safe}"

    # Rename to original plate name if sanitised name differs
    if [ "${safe}" != "${plate}" ]; then
        for ext in bed bim fam log; do
            [ -f "${safe}.\${ext}" ] && mv "${safe}.\${ext}" "${plate}.\${ext}" || true
        done
    fi
    """
}