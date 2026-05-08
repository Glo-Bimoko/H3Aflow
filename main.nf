#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    idat2vcf-pipeline  |  main.nf
    Stages:
      0.  PREP_INPUTS       – resolve idat dirs from ready/ + derive sex_info from
                              the samplesheet's "Collected Gender" column
      1.  LINK_IDATS        – symlink resolved idat directories into work dir
      2.  IDAT_TO_GTC       – bcftools +gtc2vcf plugin: idat → gtc
      3.  GTC_TO_VCF        – bcftools +gtc2vcf plugin: gtc → normalised BCF
      4.  VCF_TO_PLINK      – BCF → PLINK bed/bim/fam
      5.  GENERATE_PHENOFILE– build .phe file from samplesheet
      --- QC / sex-inference extension ---
      6.  MERGE_PLINK       – merge per-plate PLINK sets into one cohort dataset
      7.  SAMPLE_QC         – call rate + heterozygosity per sample
      8.  SPLIT_CHROM       – extract chrX and chrY beds
      9.  CHECK_SEX         – PLINK --check-sex (chrX F-statistic)
      10. XY_INTENSITY      – extract X/Y raw intensities from GTC .tsv files
      11. IBD               – PLINK --genome (IBS/IBD)
      12. PCA               – PLINK2 --pca
      13. REPORT            – compile HTML + CSV + flagged-sample list
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREP_INPUTS       } from './modules/prep_inputs'
include { LINK_IDATS        } from './modules/link_idats'
include { IDAT_TO_GTC       } from './modules/idat_to_gtc'
include { GTC_TO_VCF        } from './modules/gtc_to_vcf'
include { VCF_TO_PLINK      } from './modules/vcf_to_plink'
include { GENERATE_PHENOFILE} from './modules/generate_phenofile'
include { MERGE_PLINK       } from './modules/merge_plink'
include { SAMPLE_QC         } from './modules/sample_qc'
include { FILTER_SAMPLES   } from './modules/filter_samples'
include { SNP_QC            } from './modules/snp_qc'
include { SPLIT_CHROM       } from './modules/split_chrom'
include { CHECK_SEX         } from './modules/check_sex'
include { XY_INTENSITY      } from './modules/xy_intensity'
include { IBD               } from './modules/ibd'
include { PCA               } from './modules/pca'
include { REPORT            } from './modules/report'

// ── Validate required params ──────────────────────────────────────────────────
if (!params.samplesheet) {
    error "ERROR: --samplesheet is required.\n" +
          "  Expected columns: Sample ID, BeadChip Barcode, Sentrix Position,\n" +
          "                    Plate Number, Well Position, Collected Gender"
}

workflow {

    // -------------------------------------------------------------------------
    // Stage 0 – PREP_INPUTS
    // Reads the raw samplesheet CSV + scans ready/ recursively to resolve the
    // idat directory for each sample (matched by BeadChip Barcode + Sentrix
    // Position).  Also extracts Collected Gender → sex_info.tsv.
    // Outputs:
    //   resolved_samplesheet  – CSV with columns: sample_id, idat_dir, plate
    //   sex_info              – TSV with columns: sampleid, sex (0=F, 1=M)
    // -------------------------------------------------------------------------
    PREP_INPUTS(
        file(params.samplesheet)
    )

    // -------------------------------------------------------------------------
    // Stage 0b – parse the resolved samplesheet into a Nextflow channel
    // Each row becomes: tuple(sample_id, idat_dir_path, plate)
    // -------------------------------------------------------------------------
    ch_samplesheet = PREP_INPUTS.out.resolved_samplesheet
        .splitCsv(header: true)
        .map { row ->
            // idat_dir kept as a plain string (val) so the full absolute
            // path is preserved and not reduced to a basename by Nextflow staging
            tuple(row.sample_id, row.idat_dir, row.plate)
        }

    ch_sex_info = PREP_INPUTS.out.sex_info

    // -------------------------------------------------------------------------
    // Stage 1 – link idats
    // -------------------------------------------------------------------------
    LINK_IDATS(ch_samplesheet)

    // -------------------------------------------------------------------------
    // Stage 2 – idat → gtc  (per sample)
    // -------------------------------------------------------------------------
    IDAT_TO_GTC(
        LINK_IDATS.out.linked,
        params.bpm,   // passed as val — full path preserved
        params.egt    // passed as val — full path preserved
    )

    // -------------------------------------------------------------------------
    // Stage 3 – gtc → BCF  (per plate: group by plate)
    // -------------------------------------------------------------------------
    ch_gtcs_by_plate = IDAT_TO_GTC.out.gtc
        .map { sample_id, gtc, plate -> tuple(plate, gtc) }
        .groupTuple()

    GTC_TO_VCF(
        ch_gtcs_by_plate,
        params.bpm,
        params.egt,
        params.fasta
    )

    // -------------------------------------------------------------------------
    // Stage 4 – BCF → PLINK  (per plate)
    // -------------------------------------------------------------------------
    // Join BCF and its index before passing to VCF_TO_PLINK
    ch_bcf_with_index = GTC_TO_VCF.out.bcf
        .join(GTC_TO_VCF.out.bcf_index, by: 0)
    VCF_TO_PLINK(ch_bcf_with_index)

    // -------------------------------------------------------------------------
    // Stage 5 – phenofile  (uses the resolved samplesheet so column names
    //            are already normalised by prep_inputs.py)
    // -------------------------------------------------------------------------
    GENERATE_PHENOFILE(PREP_INPUTS.out.resolved_samplesheet)

    // -------------------------------------------------------------------------
    // Stage 6 – merge all per-plate PLINK sets → cohort dataset
    // -------------------------------------------------------------------------
    ch_plink_beds = VCF_TO_PLINK.out.plink.map { plate, bed, bim, fam -> [bed, bim, fam] }.flatten().collect()
    MERGE_PLINK(ch_plink_beds)

    // -------------------------------------------------------------------------
    // Stage 7 – sample QC: call rate + heterozygosity
    // -------------------------------------------------------------------------
    SAMPLE_QC(MERGE_PLINK.out.merged)

    // -------------------------------------------------------------------------
    // Stage 7b – SNP-level QC (geno, MAF, HWE) — mirrors h3agwas qc
    // -------------------------------------------------------------------------
    FILTER_SAMPLES(MERGE_PLINK.out.merged, SAMPLE_QC.out.keep_list)
    SNP_QC(FILTER_SAMPLES.out.plink)

    // -------------------------------------------------------------------------
    // Stage 8 – split chrX / chrY (from SNP-QC-passed data)
    // -------------------------------------------------------------------------
    SPLIT_CHROM(SNP_QC.out.plink)

    // -------------------------------------------------------------------------
    // Stage 9 – chrX F-statistic sex check
    // -------------------------------------------------------------------------
    CHECK_SEX(
        SPLIT_CHROM.out.chrX,
        ch_sex_info
    )

    // -------------------------------------------------------------------------
    // Stage 10 – X/Y raw intensity from GTC .tsv files
    // -------------------------------------------------------------------------
    ch_tsv_files = GTC_TO_VCF.out.tsv.map { plate, tsv -> tsv }.collect()
    XY_INTENSITY(ch_tsv_files, ch_sex_info)

    // -------------------------------------------------------------------------
    // Stage 11 – IBD / duplicate detection
    // -------------------------------------------------------------------------
    IBD(
        MERGE_PLINK.out.merged,
        SAMPLE_QC.out.keep_list
    )

    // -------------------------------------------------------------------------
    // Stage 12 – PCA
    // -------------------------------------------------------------------------
    PCA(
        MERGE_PLINK.out.merged,
        SAMPLE_QC.out.keep_list
    )

    // -------------------------------------------------------------------------
    // Stage 13 – report
    // -------------------------------------------------------------------------
    REPORT(
        CHECK_SEX.out.sexcheck,
        XY_INTENSITY.out.xy_tsv,
        SAMPLE_QC.out.qc_stats,
        IBD.out.genome,
        PCA.out.eigenvec,
        ch_sex_info
    )
}