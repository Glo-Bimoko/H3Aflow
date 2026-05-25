#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input  = "./data/*.txt"
params.outdir = "./results"

workflow {

    // Normalise input: if user passed a bare directory (with or without
    // trailing slash), automatically append /*.txt
    input_glob = params.input.endsWith('/') || !params.input.contains('*')
        ? "${params.input.replaceAll('/+$', '')}/*.txt"
        : params.input

    input_files = Channel.fromPath(input_glob)

    // Strip everything after the first dot: "GX099.CEL" -> "GX099"
    parsed = input_files.map { file ->
        tuple(file.getBaseName().replaceAll(/\..*/, ''), file)
    } | PARSE_SAMPLE

    // Collect all per-sample TSVs and pass them together to COMPARE.
    // No wide merge step — comparison is done pair-by-pair in compare_files.py.
    all_tsvs = parsed.map { it[1] }.collect()
    COMPARE(all_tsvs)
}

// ---------------------------------------------------------------------------
// PARSE_SAMPLE
// Extracts rsid + genotype call from a raw Axiom export .txt file.
// Output TSV has two columns: rsid (index) and <sample_id> (genotype call).
// ---------------------------------------------------------------------------
process PARSE_SAMPLE {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(file)

    output:
    tuple val(sample_id), path("${sample_id}.tsv")

    script:
    """
    python3 ${projectDir}/bin/parse_axiom.py "${file}" "${sample_id}.tsv"
    """
}

// ---------------------------------------------------------------------------
// COMPARE
// Receives all per-sample TSVs. compare_files.py loads pairs one at a time,
// joins on rsid (cross-array safe), and writes pairwise similarity results.
// No wide matrix is ever built — memory stays flat regardless of sample count.
// ---------------------------------------------------------------------------
process COMPARE {
    tag "comparing"
    publishDir params.outdir, mode: 'copy'

    input:
    path(tsvs)

    output:
    path("similarity_results.txt")

    script:
    """
    python3 ${projectDir}/bin/compare_files.py "similarity_results.txt" ${tsvs}
    """
}

