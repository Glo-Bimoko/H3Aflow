process XY_INTENSITY {
    publishDir "${params.outdir}/qc/xy_intensity", mode: 'copy'

    input:
    path tsv_files   // all per-plate .tsv files from GTC_TO_VCF
    path sex_info

    output:
    path "xy_intensity.tsv",       emit: xy_tsv
    path "xy_intensity_plot.html", emit: xy_plot

    script:
    """
    # Create header for the concatenated file
    echo -e "SAMPLE_ID\tCHR\tPOS\tREF\tALT\tNORMX\tNORMY" > all_samples.tsv
    
    # Convert Nextflow's file list to bash array and concatenate
    files=(${tsv_files})
    for f in "\${files[@]}"; do
        # Skip if it's the sex_info file or our output file
        if [[ ! "\$f" =~ sex_info ]] && [[ "\$f" != "all_samples.tsv" ]]; then
            cat "\$f" >> all_samples.tsv
        fi
    done

    python ${projectDir}/bin/extract_xy_intensity.py \
        --tsv      all_samples.tsv \
        --sex_info ${sex_info}     \
        --out      xy_intensity.tsv \
        --plot     xy_intensity_plot.html
    """
}