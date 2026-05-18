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
    # NOTE: --const-fid 0 has been removed and replaced with --double-id.
    #       make_sex_update.py writes FID=0, which would not match FID=IID.
    #       The fix is applied in check_sex.nf: the sex update file now uses
    #       FID=IID (written by make_sex_update.py after this change).
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