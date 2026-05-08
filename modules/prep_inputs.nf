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
        --out_samplesheet resolved_samplesheet.csv  \
        --out_sex_info    sex_info.tsv

    # Reset the GTC bulk-scan lock so each pipeline run gets a fresh scan
    rm -f "${params.outdir}/gtc/.scan_done"
    """

}
