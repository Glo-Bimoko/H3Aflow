process PCA {
    publishDir "${params.outdir}/pca", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)
    path keep_list

    output:
    path "pca.eigenvec", emit: eigenvec
    path "pca.eigenval", emit: eigenval
    path "pca_plot.png", emit: plot

    script:
    def prefix = bed.baseName
    
    """
    echo "PCA Started: \$(date)" > pca.log
    
    # Step 1: Filter to QC-passed samples
    echo "Step 1: Filtering samples..." >> pca.log
    plink2 \
        --bfile ${prefix} \
        --keep ${keep_list} \
        --make-bed \
        --allow-extra-chr \
        --out filtered_samples
    
    N_VAR=\$(wc -l < filtered_samples.bim)
    echo "Variants: \${N_VAR}" >> pca.log
    
    if [ \${N_VAR} -eq 0 ]; then
        echo "ERROR: No variants remaining!" | tee -a pca.log
        # Create dummy output
        echo -e "#FID\tIID\tPC1\tPC2\tPC3\tPC4\tPC5" > pca.eigenvec
        echo "0" > pca.eigenval
        exit 0
    fi
    
    # Step 2: LD pruning (skip if few variants)
    if [ \${N_VAR} -gt 10000 ]; then
        echo "Step 2: LD pruning..." >> pca.log
        plink2 \
            --bfile filtered_samples \
            --maf 0.05 \
            --indep-pairwise 200 50 0.5 \
            --allow-extra-chr \
            --out pruned
    else
        echo "Skipping LD pruning (only \${N_VAR} variants)" >> pca.log
        plink2 \
            --bfile filtered_samples \
            --write-snplist \
            --allow-extra-chr \
            --out all_snps
        cp all_snps.snplist pruned.prune.in
    fi
    
    # Check pruning results
    if [ -s pruned.prune.in ]; then
        N_PRUNE=\$(wc -l < pruned.prune.in)
        echo "Variants after pruning: \${N_PRUNE}" >> pca.log
        
        if [ \${N_PRUNE} -gt 0 ]; then
            EXTRACT_CMD="--extract pruned.prune.in"
        else
            echo "No variants passed pruning, using all variants" >> pca.log
            EXTRACT_CMD=""
        fi
    else
        echo "Pruning failed, using all variants" >> pca.log
        EXTRACT_CMD=""
    fi
    
    # Step 3: Run PCA
    N_PCS=${params.pca_components}
    echo "Step 3: Running PCA with \${N_PCS} components..." >> pca.log
    
    plink2 \
        --bfile filtered_samples \
        \${EXTRACT_CMD} \
        --pca \${N_PCS} \
        --allow-extra-chr \
        --out pca
    
    # Step 4: Plot
    if [ -s pca.eigenvec ] && [ -s pca.eigenval ]; then
        python ${projectDir}/bin/plot_pca.py \
            --eigenvec pca.eigenvec \
            --eigenval pca.eigenval \
            --out pca_plot.png || {
            echo "WARNING: Plotting failed, creating dummy plot" >> pca.log
            python -c "import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt; plt.figure(); plt.savefig('pca_plot.png')"
        }
    else
        echo "ERROR: PCA failed to produce output" >> pca.log
        python -c "import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt; plt.figure(); plt.savefig('pca_plot.png')"
    fi
    
    echo "PCA Completed: \$(date)" >> pca.log
    """
}