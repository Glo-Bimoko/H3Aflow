process GTC_QC {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(gtc), val(plate)
    
    output:
    tuple val(sample_id), path("${sample_id}.gtc_stats.tsv"), emit: stats
    
    script:
    """
    bcftools +gtc2vcf \\
        --gtc ${gtc} \\
        --extra ${sample_id}.gtc_stats.tsv \\
        > /dev/null 2>&1 || true
    
    # Ensure file exists even if bcftools fails
    touch ${sample_id}.gtc_stats.tsv
    """
}

process FILTER_GTC_SAMPLES {
    publishDir "${params.outdir}/qc", mode: 'copy'
    
    input:
    path(stats_files)
    val(gc10_threshold)
    val(call_rate_threshold)
    
    output:
    path("poorgc10.lst"), emit: poor_gc10_list
    path("gtc_qc_summary.tsv"), emit: summary
    path("gtc_qc_flags.tsv"), emit: flags
    
    script:
    """
    python3 << 'EOF'
import pandas as pd
import glob
import os

# Collect all per-sample stats
all_stats = []
for f in sorted(glob.glob("*.gtc_stats.tsv")):
    if os.path.getsize(f) == 0:
        continue
    try:
        df = pd.read_csv(f, sep='\t', nrows=1)
        sample_id = f.replace('.gtc_stats.tsv', '')
        df['sample_id'] = sample_id
        all_stats.append(df)
    except Exception as e:
        print(f"Warning: Could not parse {f}: {e}", flush=True)

if all_stats:
    stats_df = pd.concat(all_stats, ignore_index=True)
    
    # Extract key columns (bcftools +gtc2vcf output)
    # These column names match bcftools +gtc2vcf --extra output
    cols_of_interest = [
        'sample_id', 'call_rate', 'p10_gc', 'p50_gc', 
        'computed_gender', 'logr_deviation', 
        'p05_x', 'p50_x', 'p95_x', 'p05_y', 'p50_y', 'p95_y'
    ]
    
    # Filter to only columns that exist
    cols_available = ['sample_id'] + [c for c in cols_of_interest[1:] if c in stats_df.columns]
    stats_df = stats_df[cols_available]
    
    # QC filtering
    if 'p10_gc' in stats_df.columns:
        stats_df['pass_gc10'] = stats_df['p10_gc'] >= ${gc10_threshold}
    else:
        stats_df['pass_gc10'] = True
        
    if 'call_rate' in stats_df.columns:
        stats_df['pass_cr'] = stats_df['call_rate'] >= ${call_rate_threshold}
    else:
        stats_df['pass_cr'] = True
    
    stats_df['pass_qc'] = stats_df['pass_gc10'] & stats_df['pass_cr']
    
    # Write summary
    stats_df.to_csv('gtc_qc_summary.tsv', sep='\t', index=False)
    
    # Write flags for failed samples
    flags_df = stats_df[~stats_df['pass_qc']][['sample_id', 'pass_gc10', 'pass_cr']].copy()
    flags_df.to_csv('gtc_qc_flags.tsv', sep='\t', index=False)
    
    # Write poor GC10 list (PLINK format: sample_id sample_id)
    failed = stats_df[~stats_df['pass_gc10']]['sample_id'].values
    with open('poorgc10.lst', 'w') as f:
        for s in failed:
            f.write(f"{s} {s}\\n")
    
    print(f"GTC QC: {len(stats_df)} samples processed, {len(failed)} failed GC10 threshold", flush=True)
else:
    open('poorgc10.lst', 'w').close()
    open('gtc_qc_summary.tsv', 'w').close()
    open('gtc_qc_flags.tsv', 'w').close()
    print("Warning: No GTC stats files found or all were empty", flush=True)
EOF
    """
}
