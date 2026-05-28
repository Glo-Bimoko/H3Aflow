process IDAT_TO_GTC {
    tag "$sample_id"
    publishDir "${params.outdir}/gtc", mode: 'copy', pattern: "*.gtc", overwrite: true

    input:
    tuple val(sample_id), path(idat_dir), val(barcode), val(position), val(plate)
    val bpm
    val egt

    output:
    tuple val(sample_id), path("${sample_id}.gtc"), val(plate), emit: gtc

    script:
    """
    cd "\$NXF_TASK_WORKDIR"
    ${params.python} ${projectDir}/bin/convert_idat2gtc.py \
        --bpm     "${bpm}"                              \
        --egt     "${egt}"                              \
        --idats   "${idat_dir}"                         \
        --output  "\$NXF_TASK_WORKDIR/${sample_id}.gtc" \
        --gtc-dir "${params.outdir}/gtc/"
    """
}
