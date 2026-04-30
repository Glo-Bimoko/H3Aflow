process XY_INTENSITY {
    publishDir "${params.outdir}/qc/xy_intensity", mode: 'copy'

    input:
    path tsv_files   // all per-plate .tsv files from GTC_TO_VCF
    path sex_info

    output:
    path "xy_intensity.tsv",      emit: xy_tsv
    path "xy_intensity_plot.png", emit: xy_plot

    script:
    """
    # Concatenate all plate TSV files (header from first file only)
    head -1 \$(ls *.tsv | head -1) > all_samples.tsv
    for f in *.tsv; do tail -n +2 "\$f" >> all_samples.tsv; done

    python ${projectDir}/bin/extract_xy_intensity.py \
        --tsv      all_samples.tsv \
        --sex_info ${sex_info}     \
        --out      xy_intensity.tsv \
        --plot     xy_intensity_plot.png
    """
}
