process SEED_GTC {
    tag "$sample_id"
    publishDir "${params.outdir}/gtc", mode: 'copy', pattern: "*.gtc", overwrite: true

    input:
    tuple val(sample_id), val(plate)

    output:
    tuple val(sample_id), path("${sample_id}.gtc"), val(plate), emit: gtc

    script:
    """
    ${params.python} \
        ${projectDir}/bin/seed_gtc.py \
        --sample_id ${sample_id}      \
        --gtc_dir   ${params.outdir}/gtc \
        --output    ${sample_id}.gtc
    """
}