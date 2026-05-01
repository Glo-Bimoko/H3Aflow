process CHECK_SEX {
    tag "cohort"
    publishDir "${params.outdir}/qc/check_sex", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)
    path sex_info

    output:
    path "sexcheck.txt", emit: sexcheck
    path "sex_discordance.txt", emit: discordance
    path "sex_plot.png", emit: plot

    script:
    def prefix = bed.baseName
    
    """
    echo "Sex check started: \$(date)" > sex_check.log
    
    # Check if chrX file has any variants
    N_VAR=\$(wc -l < ${prefix}.bim 2>/dev/null || echo 0)
    
    if [ \${N_VAR} -eq 0 ]; then
        echo "ERROR: No chrX variants found! Cannot perform sex check." >> sex_check.log
        echo -e "IID\tCOLLECTED_SEX\tINFERRED_SEX\tF\tSTATUS\tCONCORDANT" > sexcheck.txt
        echo -e "IID\tCOLLECTED_SEX\tINFERRED_SEX\tF\tSTATUS\tCONCORDANT" > sex_discordance.txt
        python -c "import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt; plt.figure(); plt.savefig('sex_plot.png')"
        exit 0
    fi
    
    # Step 1: Run PLINK's sex check on chrX
    plink2 \
        --bfile ${prefix} \
        --check-sex \
        --allow-extra-chr \
        --out sexcheck 2>&1 | tee -a sex_check.log
    
    # Step 2: Parse results and compare with collected sex
    if [ -f sexcheck.sexcheck ]; then
        python ${projectDir}/bin/compare_sex.py \
            --plink_sex sexcheck.sexcheck \
            --collected_sex ${sex_info} \
            --out sex_discordance.txt \
            --plot sex_plot.png 2>&1 | tee -a sex_check.log
        
        # Create sexcheck.txt as a simpler version of the output
        if [ -f sex_discordance.txt ]; then
            cp sex_discordance.txt sexcheck.txt
        else
            echo -e "IID\tCOLLECTED_SEX\tINFERRED_SEX\tF\tSTATUS\tCONCORDANT" > sexcheck.txt
        fi
    else
        echo "ERROR: PLINK sex check failed to produce output" >> sex_check.log
        echo -e "IID\tCOLLECTED_SEX\tINFERRED_SEX\tF\tSTATUS\tCONCORDANT" > sexcheck.txt
        echo -e "IID\tCOLLECTED_SEX\tINFERRED_SEX\tF\tSTATUS\tCONCORDANT" > sex_discordance.txt
        python -c "import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt; plt.figure(); plt.savefig('sex_plot.png')"
    fi
    
    echo "Sex check completed: \$(date)" >> sex_check.log
    """
}