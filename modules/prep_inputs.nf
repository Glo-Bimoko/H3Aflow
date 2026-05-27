process PREP_INPUTS {
    publishDir "${params.outdir}/inputs", mode: 'copy'

    input:
    path samplesheet   // raw project samplesheet CSV

    output:
    path "resolved_samplesheet.csv", emit: resolved_samplesheet
    path "sex_info.tsv",             emit: sex_info

    script:
    """
    python ${projectDir}/bin/prep_inputs.py \
        --samplesheet     ${samplesheet}            \
        --idat_root       ${params.idat_root}       \
        --gtc_dir         ${params.outdir}/gtc      \
        --out_samplesheet resolved_samplesheet.csv  \
        --out_sex_info    sex_info.tsv              \
        --reset-gtc-scan  --prefer-existing-gtc
    """

}
