process IDAT_TO_GTC {
    tag "$sample_id"
    publishDir "${params.outdir}/gtc", mode: 'copy', pattern: "*.gtc"

    input:
    tuple val(sample_id), path(idat_dir), val(plate)
    val bpm
    val egt

    output:
    tuple val(sample_id), path("${sample_id}.gtc"), val(plate), emit: gtc

    script:
    """
    python ${projectDir}/bin/convert_idat2gtc.py \
        --bpm    "${bpm}"       \
        --egt    "${egt}"       \
        --idats  "${idat_dir}"  \
        --output ${sample_id}.gtc
    """
}

