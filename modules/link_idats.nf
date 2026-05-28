process LINK_IDATS {
    tag "$sample_id"
    publishDir "${params.outdir}/linked_idats", mode: 'copy'

    input:
    tuple val(sample_id), val(idat_dir), val(barcode), val(position), val(plate)

    output:
    tuple val(sample_id), path("linked_idats/${sample_id}"), val(barcode), val(position), val(plate), emit: linked

    script:
    """
    ${params.python} \
        ${projectDir}/bin/link_idats.py \
        --source   "${idat_dir}"        \
        --dest     linked_idats/${sample_id} \
        --barcode  "${barcode}"         \
        --position "${position}"
    """
}