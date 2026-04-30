process LINK_IDATS {
    tag "$sample_id"

    input:
    tuple val(sample_id), val(idat_dir), val(plate)

    output:
    tuple val(sample_id), path("linked_idats/${sample_id}"), val(plate), emit: linked

    script:
    """
    python ${projectDir}/bin/link_idats.py \
        --source "${idat_dir}" \
        --dest   linked_idats/${sample_id}
    """
}

