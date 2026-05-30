# H3Aflow

**A Nextflow DSL2 pipeline for end-to-end processing of Illumina H3Africa SNP array data. From raw idat or GTC files to QC-passed, analysis-ready PLINK datasets.**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10-23aa62.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Overview

H3Aflow bridges the gap between raw Illumina array files and the PLINK-format input required by downstream GWAS tools such as [h3agwas](https://github.com/h3abionet/h3agwas). No existing public pipeline covered this full conversion path for H3Africa array data. H3Aflow was built to fill that gap.

```
idat files  ──OR──  GTC files (pre-existing)
    │                    │
    │                    └─► SEED_GTC (copy from results/gtc/ → skip conversion)
    │
    ▼
Stage 0   PREP_INPUTS        Resolve idat dirs from samplesheet · detect existing GTCs · derive sex_info
Stage 1   LINK_IDATS         Symlink per-sample idat files into clean work dirs
Stage 2   IDAT_TO_GTC        bcftools +idat2gtc  →  GTC files (per sample, skipped if GTC provided)
Stage 2a  GTC_QC             bcftools +gtc2vcf stats  →  per-plate GTC quality metrics
Stage 3   GTC_TO_VCF         bcftools +gtc2vcf   →  normalised BCF + XY intensity TSV (per plate)
Stage 4   VCF_TO_PLINK       plink --bcf         →  PLINK bed/bim/fam (per plate)
Stage 5   GENERATE_PHENOFILE Build PLINK .phe file from samplesheet
Stage 6   MERGE_PLINK        Merge all per-plate PLINK sets → one cohort dataset
Stage 7   SAMPLE_QC          Flag samples: call rate (--mind) + heterozygosity outliers
Stage 7b  SNP_QC             Filter SNPs: missingness (--geno) · MAF · HWE (autosomes only)
Stage 8   CHECK_SEX          GTC computed_gender vs collected sex + plate discordance report
Stage 9   GENOTYPE_CONCORDANCE  All-vs-all pairwise concordance · flag likely duplicates
Stage 10  IBD                LD pruning → plink --genome → flag duplicate/related pairs
Stage 11  PCA                plink2 --pca → population structure
Stage 12  REPORT             Self-contained HTML QC report + flagged_samples.tsv
    │
    ▼
cohort.bed / .bim / .fam   ←  ready for h3agwas or association testing
```

---

## Required Starting Materials

H3Aflow accepts either **idat files** or **GTC files** as input. One does not need both.

**Option A — idat files (raw fluorescence)**
Place all idat files under a single root directory (they can be nested arbitrarily deep). The pipeline discovers them by matching `<BeadChip Barcode>_<Sentrix Position>_Red.idat` and `_Grn.idat` patterns recursively and converts them to GTC format automatically.

**Option B — GTC files (pre-converted)**
If you already have GTC files from a previous conversion run, place them in `results/gtc/` using the naming convention `{sample_id}.gtc`. The pipeline detects them automatically during `PREP_INPUTS`, routes those samples through `SEED_GTC` (a fast copy step), and skips the `IDAT_TO_GTC` conversion entirely for those samples.

Both input types can be mixed. Samples with existing GTCs skip conversion while samples with only idat files go through the full conversion path. This is particularly useful when resuming a partially completed run or when integrating data from multiple sources.

---

## Requirements

| Tool | Version | Purpose |
|---|---|---|
| [Nextflow](https://www.nextflow.io/) | ≥ 22.10 | Workflow engine |
| [bcftools](https://samtools.github.io/bcftools/) | ≥ 1.11 | idat→GTC, GTC→VCF |
| [bcftools gtc2vcf plugin](https://github.com/freeseek/gtc2vcf) | 2024-05-05 | idat2gtc, gtc2vcf .so files |
| [PLINK 1.9](https://www.cog-genomics.org/plink/) | 1.9 | VCF→PLINK, sample QC, IBD |
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
git clone https://github.com/Glo-Bimoko/H3Aflow.git
cd H3Aflow

# Install Python dependencies into your conda environment
pip install pandas numpy matplotlib openpyxl

# Set bcftools plugin path
export BCFTOOLS_PLUGINS=/path/to/gtc2vcf/plugins
```

---

## Quick Start

### 1. Prepare your samplesheet

The samplesheet is a CSV file with the following columns (header names are case-insensitive, extra columns are ignored):

| Column | Description | Example |
|---|---|---|
| Sample ID | Unique sample identifier | `7801848` |
| BeadChip Barcode | Illumina array barcode used to locate idat files | `205695340002` |
| Sentrix Position | Position on the BeadChip | `R01C01` |
| Plate Number | Batch/plate grouping | `Plate_01` |
| Well Position | Optional, carried through for traceability | `A01` |
| Collected Gender | Reported sex: Female/Male, F/M, 0/1 | `Female` |

Sex information is extracted automatically from the `Collected Gender` column. No separate sex_info file is needed.

### 2. Organise your input files

**idat files:**
```
ready/
├── Plate_01_IDAT_Files/
│   ├── 205695340002_R01C01_Red.idat
│   ├── 205695340002_R01C01_Grn.idat
│   └── ...
└── Plate_02_IDAT_Files/
    └── ...
```

**GTC files (if skipping idat conversion):**
```
results/
└── gtc/
    ├── 7801848.gtc
    ├── 7801849.gtc
    └── ...
```

GTC files must be named `{sample_id}.gtc` where `sample_id` matches the Sample ID column in your samplesheet. The pipeline detects them automatically — no extra flags required.

### 3. Configure

Edit `nextflow.config` to set your paths. The `python` parameter must point to the Python interpreter in your conda environment:

```groovy
// local profile
params {
    python    = "/home/<user>/miniforge3/envs/<env_name>/bin/python"
    bpm       = "/path/to/H3Africa_2019_20037295_A1.bpm"
    egt       = "/path/to/H3Africa_2019_Gentrain_A1_ClusterFile_Final.egt"
    fasta     = "/path/to/human_g1k_v37.fasta"
    idat_root = "/path/to/ready/"
    outdir    = "/path/to/results"
}
```

### 4. Run

**Locally (desktop/laptop):**
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
| `local` | Laptop or desktop — no module system, uses local conda env |
| `chpc` | CHPC Lengau HPC — loads modules via `module load`, uses Lustre conda env |
| `pbs` | PBS GPU queue (alternative CHPC submission) |

Each profile sets its own `params.python` pointing to the correct conda environment for that machine. This means the same pipeline code runs on your desktop and on the cluster without any path changes in your `.nf` files.

### CHPC Lustre conda env auto-sync

Lustre storage at CHPC is purged every 90 days. Seeing that the scheduler isn't able to access a user's home directory, all conda envs belonging to the user will not be visible to the pipeline when a job is submitted. The `run_pipeline.qsub` script handles this automatically by making a clone of the conda environment. In the case where it was purged after 90 days, the .qsub script will auto-sync, thus removing the need for regular maintenance. At the start of each job it checks whether the conda env exists on Lustre and clones it from the user's home directory if previously purged. No manual intervention is needed.

---

## Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `--samplesheet` | required | Path to samplesheet CSV |
| `--idat_root` | required | Root directory containing idat files |
| `--bpm` | required | Illumina BPM manifest |
| `--egt` | required | Illumina EGT cluster file |
| `--fasta` | required | Reference genome FASTA |
| `--python` | `"python"` | Path to Python interpreter (set per profile in config) |
| `--mind` | `0.05` | Sample missingness cutoff (>5% → fail) |
| `--het_sd` | `3` | Heterozygosity SD cutoff |
| `--snp_missing` | `0.05` | SNP missingness cutoff (`--geno`) |
| `--maf` | `0.0001` | Minor allele frequency filter |
| `--hwe` | `1e-4` | Hardy-Weinberg equilibrium p-value cutoff |
| `--gc10_threshold` | `0.15` | Minimum GTC p10 GenCall score |
| `--call_rate_threshold` | `0.95` | Minimum GTC-level call rate |
| `--ibd_pi_hat` | `0.1875` | PI_HAT threshold for flagging related pairs |
| `--pca_components` | `10` | Number of PCs to compute |
| `--outdir` | `./results` | Output directory |

---

## Output

```
results/
├── inputs/
│   ├── resolved_samplesheet.csv   # idat paths resolved + input_source per sample
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
│   ├── sample_qc/                 # call rate + het stats + pass/fail lists
│   ├── snp_qc/                    # SNP filtering logs (autosomes only)
│   ├── check_sex/                 # GTC computed_gender vs collected sex
│   ├── concordance/               # pairwise genotype concordance
│   ├── ibd/                       # related pair flags
│   └── pca/                       # eigenvectors + scree plot
└── report/
    ├── qc_report.html             # self-contained HTML QC report
    ├── flagged_samples.tsv        # union of all QC failures
    └── cohort_summary.tsv         # cohort-level summary statistics
```

---

## Pipeline Structure

```
H3Aflow/
├── main.nf                        # Orchestrates the full workflow
├── nextflow.config                # Parameters, resource limits, profiles
├── run_pipeline.qsub              # PBS job script for CHPC Lengau
├── README.md
├── modules/
│   ├── prep_inputs.nf
│   ├── link_idats.nf
│   ├── idat_to_gtc.nf
│   ├── seed_gtc.nf
│   ├── gtc_qc.nf
│   ├── gtc_to_vcf.nf
│   ├── vcf_to_plink.nf
│   ├── generate_phenofile.nf
│   ├── merge_plink.nf
│   ├── sample_qc.nf
│   ├── filter_samples.nf
│   ├── snp_qc.nf
│   ├── check_sex.nf
│   ├── genotype_concordance.nf
│   ├── ibd.nf
│   ├── pca.nf
│   └── report.nf
└── bin/
    ├── prep_inputs.py
    ├── link_idats.py
    ├── seed_gtc.py
    ├── convert_idat2gtc.py
    ├── generate_phenofile.py
    ├── compute_sample_qc.py
    ├── pairwise_concordance.py
    ├── check_idat_duplicates.py
    ├── gtc_sex_check.py
    ├── flag_ibd_duplicates.py
    ├── plot_pca.py
    └── generate_report.py
```

---

## Key Design Decisions

**Sex check uses GTC computed_gender, not PLINK --check-sex.**
Illumina's genotype calling assigns a computed gender to each sample during GTC generation. H3Aflow extracts this from `bcftools +gtc2vcf` stats and compares it to the samplesheet's Collected Gender. This avoids the need for chrX LD pruning and is more robust on arrays with sparse chrX coverage.

**SNP QC is autosomes-only.**
`--geno`, `--maf`, and `--hwe` are applied with `--not-chr X Y XY`. Sex chromosome SNPs are excluded from these filters because male hemizygosity inflates chrX missingness and makes cohort-wide HWE invalid. chrX QC is handled separately in CHECK_SEX.

**IBD and PCA use the pre-SNP-QC merged dataset.**
Running IBD and PCA on the full variant set (before SNP QC) gives more stable estimates. The SAMPLE_QC keep-list is applied to restrict to QC-passing samples.

---

## Comparison with H3AGWAS

H3Aflow is a complete, standalone pipeline. It does not require H3AGWAS to perform QC. The two pipelines are complementary. H3Aflow covers everything from raw array files through QC, while H3AGWAS covers association testing and meta-analysis:

| Capability | H3Aflow | H3AGWAS |
|---|---|---|
| idat → GTC → VCF conversion | ✅ | ❌ |
| GTC file input (skip conversion) | ✅ | ❌ |
| H3Africa BPM/EGT support | ✅ | ❌ |
| XY raw intensity QC | ✅ | ❌ |
| Sex inference (GTC computed gender) | ✅ | Partial |
| GTC-level call rate / GenCall QC | ✅ | ❌ |
| Sample QC (call rate, het) | ✅ | ✅ |
| SNP QC (MAF, HWE, geno) | ✅ | ✅ |
| Genotype concordance | ✅ | ❌ |
| IBD / duplicate detection | ✅ | ✅ |
| PCA | ✅ | ✅ |
| Association testing | ❌ | ✅ |
| Meta-analysis | ❌ | ✅ |
| DSL2 | ✅ | ❌ (requires NXF ≤ 22.10) |

H3Aflow produces a fully QC-passed `cohort.bed / .bim / .fam` dataset that is directly compatible with the h3agwas `assoc` workflow — no format conversion needed.

---

## Acknowledgements

- [bcftools gtc2vcf plugin](https://github.com/freeseek/gtc2vcf) by Giulio Genovese (freeseek)
- [h3agwas](https://github.com/h3abionet/h3agwas) by H3ABioNet for QC design patterns
- Dr Ayoub Ksouri (https://za.linkedin.com/in/ayoub-ksouri) for providing conversion scripts used during development
- Dr. Jean-Tristan Brandenburg (https://za.linkedin.com/in/brandenburgj) for ideas and inspiration drawn from the H3AGWAS pipeline
- [CHPC](https://www.chpc.ac.za/) Lengau cluster for compute resources

---

## Citation

If you use H3Aflow in your research, please cite:

> Glory Bimoko & the CPGR Team. H3Aflow: end-to-end Nextflow pipeline for Illumina H3Africa SNP array processing. https://github.com/Glo-Bimoko/H3Aflow

---

## License

MIT License — see [LICENSE](LICENSE) for details.