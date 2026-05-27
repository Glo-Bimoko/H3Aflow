process GENERATE_PHENOFILE {
    publishDir "${params.outdir}/phenofile", mode: 'copy'

    input:
    path samplesheet

    output:
    path "sample.phe", emit: phenofile

    script:
    """
    python ${projectDir}/bin/generate_phenofile.py \
        --samplesheet ${samplesheet} \
        --out sample.phe
    """
}
