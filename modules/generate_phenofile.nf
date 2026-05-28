process GENERATE_PHENOFILE {
    publishDir "${params.outdir}/phenofile", mode: 'copy'

    input:
    path samplesheet

    output:
    path "sample.phe", emit: phenofile

    script:
    """
    ${params.python} \
        ${projectDir}/bin/generate_phenofile.py \
        --samplesheet ${samplesheet}            \
        --out sample.phe
    """
}