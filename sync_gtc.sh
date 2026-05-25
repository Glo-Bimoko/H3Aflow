#!/bin/bash

# Path to your samplesheet
SAMPLESHEET="/home/eiko/Desktop/H3Aflow/GRCh37/made_samplesheet.csv"

# Build rsync filter rules from Sample IDs (skip header line)
FILTERS=()
while IFS=',' read -r sample_id _; do
    FILTERS+=(--include="${sample_id}.gtc")
done < <(tail -n +2 "$SAMPLESHEET")

# Exclude everything that didn't match
FILTERS+=(--exclude='*')

rsync -avzP \
  "${FILTERS[@]}" \
  -e "ssh -o ServerAliveInterval=60 -o ServerAliveCountMax=5" \
  gbimoko@scp.chpc.ac.za:/home/gbimoko/lustre/idat2vcf-pipeline/idat2vcf-pipeline/results/gtc/ \
  /home/eiko/Desktop/H3Aflow/results/gtc/
