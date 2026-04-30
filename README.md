# h3aflow

**A Nextflow DSL2 pipeline for end-to-end processing of Illumina H3Africa SNP array data — from raw idat files to QC-passed, analysis-ready PLINK datasets.**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10-23aa62.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Overview

h3aflow bridges the gap between raw Illumina idat files and the PLINK-format input required by downstream GWAS tools such as [h3agwas](https://github.com/h3abionet/h3agwas). No existing public pipeline covered this full conversion path for H3Africa array data. H3aflow was built to fill that gap.

```
idat files (raw fluorescence)
    │
    ▼
Stage 0  PREP_INPUTS        Resolve idat dirs from samplesheet · derive sex_info
Stage 1  LINK_IDATS         Symlink per-sample idat files into clean work dirs
Stage 2  IDAT_TO_GTC        bcftools +idat2gtc  →  GTC files (per sample)
Stage 3  GTC_TO_VCF         bcftools +gtc2vcf   →  normalised BCF + XY intensity TSV (per plate)
Stage 4  VCF_TO_PLINK       plink --bcf         →  PLINK bed/bim/fam (per plate)
Stage 5  GENERATE_PHENOFILE Build PLINK .phe file from samplesheet
Stage 6  MERGE_PLINK        Merge all per-plate PLINK sets → one cohort dataset
Stage 7  SAMPLE_QC          Flag samples: call rate (--mind) + heterozygosity outliers
Stage 7b SNP_QC             Filter SNPs: missingness (--geno) · MAF · HWE
Stage 8  SPLIT_CHROM        Extract chrX and chrY for sex inference
Stage 9  CHECK_SEX          plink --check-sex (chrX F-statistic) + collected sex comparison
Stage 10 XY_INTENSITY       Extract raw X/Y intensities from GTC TSVs
Stage 11 IBD                LD pruning → plink --genome → flag duplicate/related pairs
Stage 12 PCA                plink2 --pca → population structure
Stage 13 REPORT             Self-contained HTML QC report + flagged_samples.tsv
    │
    ▼
cohort_qc.bed / .bim / .fam   ←  ready for h3agwas qc or association testing
```

---

## Requirements

| Tool | Version | Purpose |
|---|---|---|
| [Nextflow](https://www.nextflow.io/) | ≥ 22.10 | Workflow engine |
| [bcftools](https://samtools.github.io/bcftools/) | ≥ 1.11 | idat→GTC, GTC→VCF |
| [bcftools gtc2vcf plugin](https://github.com/freeseek/gtc2vcf) | 2024-05-05 | idat2gtc, gtc2vcf .so files |
| [PLINK 1.9](https://www.cog-genomics.org/plink/) | 1.9 | VCF→PLINK, sample QC, sex check, IBD |
| [PLINK 2](https://www.cog-genomics.org/plink/2.0/) | ≥ 2.0 | PCA |
| Python | ≥ 3.8 | QC scripts |
| pandas, numpy, matplotlib | latest | Python dependencies |

### Array assets required

You must supply your own Illumina manifest and cluster files matching your array version:

| File | Description |
|---|---|
| `*.bpm` | Illumina BPM manifest file |
| `*.egt` | Illumina EGT cluster file |
| `*.fasta` + `.fai` | Reference genome (GRCh37 recommended for H3Africa arrays) |

> **Important:** The BPM and EGT must be from the same array revision. The H3Africa 2019 array has A1 and B1 revisions — mismatching them will cause probe lookup failures.

---

## Installation

```bash
git clone https://github.com/gbimoko/h3aflow.git
cd h3aflow

# Install Python dependencies
pip install pandas numpy matplotlib openpyxl

# Set bcftools plugin path
export BCFTOOLS_PLUGINS=/path/to/gtc2vcf/plugins
```

---

## Quick Start

### 1. Prepare your samplesheet

The samplesheet is a CSV file with the following columns (header names are case-insensitive):

| Column | Description | Example |
|---|---|---|
| Sample ID | Unique sample identifier | `7801848` |
| BeadChip Barcode | Illumina array barcode — used to locate idat files | `205695340002` |
| Sentrix Position | Position on the BeadChip | `R01C01` |
| Plate Number | Batch/plate grouping | `Plate_01` |
| Well Position | Optional, carried through for traceability | `A01` |
| Collected Gender | Reported sex: Female/Male, F/M, 0/1 | `Female` |

> Sex information is extracted automatically from the `Collected Gender` column. No separate sex_info file is needed.

### 2. Organise your idat files

Place all idat files under a single root directory (they can be nested arbitrarily deep):

```
ready/
├── Plate_01_IDAT_Files/
│   ├── 205695340002_R01C01_Red.idat
│   ├── 205695340002_R01C01_Grn.idat
│   └── ...
└── Plate_02_IDAT_Files/
    └── ...
```

The pipeline discovers idat files by matching `<BeadChip Barcode>_<Sentrix Position>_Red.idat` and `_Grn.idat` patterns recursively.

### 3. Configure

Edit `nextflow.config` to set your paths:

```groovy
params {
    bpm      = "/path/to/H3Africa_2019_20037295_A1.bpm"
    egt      = "/path/to/H3Africa_2019_Gentrain_A1_ClusterFile_Final.egt"
    fasta    = "/path/to/human_g1k_v37.fasta"
    idat_root = "/path/to/ready/"
    outdir   = "/path/to/results"
}
```

### 4. Run

**Locally:**
```bash
nextflow run main.nf \
  --samplesheet /path/to/samplesheet.csv \
  -profile local \
  -resume
```

**On CHPC Lengau (PBS):**
```bash
qsub run_pipeline.qsub
```

---

## Configuration Profiles

| Profile | Use case |
|---|---|
| `local` | Laptop or desktop, no module system |
| `chpc` | CHPC Lengau HPC, loads modules via `module load` |
| `pbs` | PBS job submission (future use) |

---

## Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `--samplesheet` | required | Path to samplesheet CSV |
| `--idat_root` | required | Root directory containing idat files |
| `--bpm` | required | Illumina BPM manifest |
| `--egt` | required | Illumina EGT cluster file |
| `--fasta` | required | Reference genome FASTA |
| `--mind` | `0.05` | Sample missingness cutoff (>5% → fail) |
| `--het_sd` | `3` | Heterozygosity SD cutoff |
| `--snp_missing` | `0.05` | SNP missingness cutoff (`--geno`) |
| `--maf` | `0.01` | Minor allele frequency filter |
| `--hwe` | `0.00001` | Hardy-Weinberg equilibrium p-value cutoff |
| `--ibd_pi_hat` | `0.1875` | PI_HAT threshold for flagging related pairs |
| `--pca_components` | `10` | Number of PCs to compute |
| `--outdir` | `./results` | Output directory |

---

## Output

```
results/
├── inputs/
│   ├── resolved_samplesheet.csv   # idat paths resolved for each sample
│   └── sex_info.tsv               # sex coded 0=Female, 1=Male
├── gtc/                           # per-sample GTC files
├── bcf/                           # per-plate normalised BCF files
├── tsv/                           # per-plate raw XY intensity tables
├── plink_per_plate/               # per-plate PLINK bed/bim/fam
├── plink_merged/
│   └── cohort.bed / .bim / .fam  # full cohort PLINK dataset
├── phenofile/
│   └── sample.phe                 # PLINK phenotype file
├── qc/
│   ├── sample_qc/                 # call rate + het stats
│   ├── snp_qc/                    # SNP filtering logs
│   ├── sex_check/                 # F-statistic + discordant samples
│   ├── xy_intensity/              # XY intensity plots
│   ├── ibd/                       # related pair flags
│   └── pca/                       # eigenvectors + scree plot
└── report/
    ├── qc_report.html             # self-contained HTML QC report
    ├── flagged_samples.tsv        # union of all QC failures
    └── cohort_summary.tsv         # cohort-level summary statistics
```

---

## Comparison with h3agwas

h3aflow is a complete, standalone pipeline. It does not require h3agwas to perform QC. The two pipelines are complementary — h3aflow covers everything from raw idat files through QC, while h3agwas covers association testing and meta-analysis:

| Capability | h3aflow | h3agwas |
|---|---|---|
| idat → GTC → VCF conversion | ✅ | ❌ |
| H3Africa BPM/EGT support | ✅ | ❌ |
| XY raw intensity QC | ✅ | ❌ |
| Sex inference (F-stat + XY) | ✅ | Partial |
| Sample QC (call rate, het) | ✅ | ✅ |
| SNP QC (MAF, HWE, geno) | ✅ | ✅ |
| IBD / duplicate detection | ✅ | ✅ |
| PCA | ✅ | ✅ |
| Association testing | ❌ | ✅ |
| Meta-analysis | ❌ | ✅ |
| DSL2 | ✅ | ❌ (requires NXF ≤ 22.10) |

h3aflow produces a fully QC-passed `cohort_qc.bed / .bim / .fam` dataset. If you wish to proceed to association testing, this output is directly compatible with the h3agwas `assoc` workflow — no format conversion needed.

---

## Pipeline Structure

```
h3aflow/
├── main.nf                   # Orchestrates all 13 stages
├── nextflow.config           # Parameters, resource limits, profiles
├── run_pipeline.qsub         # PBS job script for CHPC Lengau
├── README.md
├── modules/
│   ├── prep_inputs.nf
│   ├── link_idats.nf
│   ├── idat_to_gtc.nf
│   ├── gtc_to_vcf.nf
│   ├── vcf_to_plink.nf
│   ├── generate_phenofile.nf
│   ├── merge_plink.nf
│   ├── sample_qc.nf
│   ├── snp_qc.nf
│   ├── split_chrom.nf
│   ├── check_sex.nf
│   ├── xy_intensity.nf
│   ├── ibd.nf
│   ├── pca.nf
│   └── report.nf
└── bin/
    ├── prep_inputs.py
    ├── link_idats.py
    ├── convert_idat2gtc.py
    ├── convert_gtc2vcf.py
    ├── generate_phenofile.py
    ├── compute_sample_qc.py
    ├── annotate_sex_check.py
    ├── extract_xy_intensity.py
    ├── flag_ibd_duplicates.py
    ├── plot_pca.py
    └── generate_report.py
```

---

## Citation

If you use h3aflow in your research, give props to:

> Glory Bimoko & the CPGR Team for making this accessible.

---

## Acknowledgements

- [bcftools gtc2vcf plugin](https://github.com/freeseek/gtc2vcf) by Giulio Genovese (freeseek)
- [h3agwas](https://github.com/h3abionet/h3agwas) by H3ABioNet for QC design patterns
- [CHPC](https://www.chpc.ac.za/) Lengau cluster for compute resources
- The H3Africa consortium for array design and data

---

## License

MIT License — see [LICENSE](LICENSE) for details.
