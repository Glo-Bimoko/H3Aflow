#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    H3Aflow  |  main.nf

    Stages:
      0.  PREP_INPUTS           – resolve idat dirs from ready/ + derive sex_info from
                                  the samplesheet's "Collected Gender" column
      1.  LINK_IDATS            – symlink resolved idat directories into work dir
      2.  IDAT_TO_GTC           – bcftools +gtc2vcf plugin: idat → gtc
      3.  GTC_TO_VCF            – bcftools +gtc2vcf plugin: gtc → normalised BCF
      4.  VCF_TO_PLINK          – BCF → PLINK bed/bim/fam
      5.  GENERATE_PHENOFILE    – build .phe file from samplesheet
      6.  MERGE_PLINK           – merge per-plate PLINK sets into one cohort dataset
      7.  SAMPLE_QC             – call rate + heterozygosity per sample
      8.  FILTER_SAMPLES        – remove failing samples

      ── Sex check (GTC computed_gender vs collected sex) ─────────────────────────
      9.  CHECK_SEX             – compare Illumina GTC computed_gender with samplesheet
                                  + per-plate discordance report + XY intensity summary

      ── SNP QC branch (autosomes only) ───────────────────────────────────────────
      11. SNP_QC                – geno / MAF / HWE on autosomes only
                                  (chrX and chrY excluded — see snp_qc.nf)

      ── Downstream QC (all depend on SNP-QC-passed autosomal dataset) ───────────
      12. GENOTYPE_CONCORDANCE  – all-vs-all pairwise genotype concordance (chunked)
      13. IBD                   – PLINK --genome (IBS/IBD)
      14. PCA                   – PLINK2 --pca
      15. REPORT                – compile HTML + CSV + flagged-sample list

    DAG overview:
      MERGE ─► SAMPLE_QC ─► FILTER ──┬─ SNP_QC ─► CONCORDANCE ─┐
                                     │                          ├─► REPORT
      GTC_QC ─► FILTER_GTC ─► CHECK_SEX (computed_gender) ─────┘
      (IBD and PCA run from MERGE+SAMPLE_QC output directly, not from SNP_QC)

    Key design decisions:
      1. Sex check uses Illumina GTC computed_gender (from bcftools +gtc2vcf stats
         during GTC_QC) compared to Collected Gender on the samplesheet — no chrX
         PLINK --check-sex branch.

      2. SNP_QC applies --not-chr X Y XY so only autosomes are filtered.

      3. CHECK_SEX emits sexcheck_plate_report.tsv for batch-level label-swap alerts.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREP_INPUTS           } from './modules/prep_inputs'
include { LINK_IDATS            } from './modules/link_idats'
include { IDAT_TO_GTC           } from './modules/idat_to_gtc'
include { GTC_QC; FILTER_GTC_SAMPLES } from './modules/gtc_qc'
include { SEED_GTC              } from './modules/seed_gtc'
include { GTC_TO_VCF            } from './modules/gtc_to_vcf'
include { VCF_TO_PLINK          } from './modules/vcf_to_plink'
include { GENERATE_PHENOFILE    } from './modules/generate_phenofile'
include { MERGE_PLINK           } from './modules/merge_plink'
include { SAMPLE_QC             } from './modules/sample_qc'
include { FILTER_SAMPLES        } from './modules/filter_samples'
include { CHECK_SEX             } from './modules/check_sex'
include { SNP_QC                } from './modules/snp_qc'
include { GENOTYPE_CONCORDANCE  } from './modules/genotype_concordance'
include { IBD                   } from './modules/ibd'
include { PCA                   } from './modules/pca'
include { REPORT                } from './modules/report'

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
    // Stage 0b – parse the resolved samplesheet into a per-sample channel
    // Each row becomes: tuple(sample_id, idat_dir, barcode, position, plate)
    // idat_dir is kept as val (not file) so the full absolute path is preserved
    // and not reduced to a basename by Nextflow staging.
    // -------------------------------------------------------------------------
    ch_samplesheet = PREP_INPUTS.out.resolved_samplesheet
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample_id,
                row.idat_dir,
                row.barcode,
                row.position,
                row.plate,
                row.input_source ?: "idat"
            )
        }

    ch_sex_info = PREP_INPUTS.out.sex_info

    ch_idat_samples = ch_samplesheet
        .filter { _sid, _idir, _bc, _pos, _plate, source -> source == "idat" }
        .map { sid, idir, bc, pos, plate, _source -> tuple(sid, idir, bc, pos, plate) }

    ch_gtc_seed_samples = ch_samplesheet
        .filter { _sid, _idir, _bc, _pos, plate, source -> source == "gtc" }
        .map { sid, _idir, _bc, _pos, plate, _source -> tuple(sid, plate) }

    // -------------------------------------------------------------------------
    // Stage 1 – LINK_IDATS
    // Symlinks each sample's idat directory into the Nextflow work directory
    // so downstream processes have stable paths.
    // -------------------------------------------------------------------------
    LINK_IDATS(ch_idat_samples)

    // -------------------------------------------------------------------------
    // Stage 2 – IDAT_TO_GTC  (per sample)
    // Calls bcftools +gtc2vcf to convert each idat pair → GTC file.
    // -------------------------------------------------------------------------
    IDAT_TO_GTC(
        LINK_IDATS.out.linked,
        params.bpm,   // val — full path preserved
        params.egt    // val — full path preserved
    )

    // -------------------------------------------------------------------------
    // Stage 2b – SEED_GTC  (pre-supplied GTC, no idat)
    // Copies {sample_id}.gtc from results/gtc when prep_inputs routes input_source=gtc.
    // -------------------------------------------------------------------------
    SEED_GTC(ch_gtc_seed_samples)

    // All GTCs entering VCF conversion (from idat conversion and/or seeding).
    ch_all_gtcs = IDAT_TO_GTC.out.gtc.mix(SEED_GTC.out.gtc)

    // -------------------------------------------------------------------------
    // Stage 2a – GTC_QC  (per plate)
    // Runs bcftools +gtc2vcf on all GTCs for a plate (same batching as GTC_TO_VCF).
    // Must use ch_all_gtcs — when every sample is SEED_GTC-only, IDAT_TO_GTC is empty.
    // -------------------------------------------------------------------------
    ch_gtcs_for_qc = ch_all_gtcs
        .map { sample_id, gtc, plate -> tuple(plate, gtc) }
        .groupTuple()

    GTC_QC(
        ch_gtcs_for_qc,
        params.bpm,
        params.egt,
        params.fasta
    )

    FILTER_GTC_SAMPLES(
        GTC_QC.out.stats.collect(),
        params.gc10_threshold ?: 0.15,
        params.call_rate_threshold ?: 0.95
    )

    // -------------------------------------------------------------------------
    // Stage 3 – GTC_TO_VCF  (per plate)
    // Groups GTC files by plate, then calls bcftools +gtc2vcf to produce a
    // normalised, reference-aligned BCF per plate.
    // -------------------------------------------------------------------------
    ch_gtcs_by_plate = ch_all_gtcs
        .map { sample_id, gtc, plate -> tuple(plate, gtc) }
        .groupTuple()

    GTC_TO_VCF(
        ch_gtcs_by_plate,
        params.bpm,
        params.egt,
        params.fasta,
        PREP_INPUTS.out.sex_info
    )

    // -------------------------------------------------------------------------
    // Stage 4 – VCF_TO_PLINK  (per plate)
    // Converts each per-plate BCF to PLINK bed/bim/fam.
    // BCF index joined before passing so both files are staged together.
    // -------------------------------------------------------------------------
    ch_bcf_with_index = GTC_TO_VCF.out.bcf
        .join(GTC_TO_VCF.out.bcf_index, by: 0)

    VCF_TO_PLINK(ch_bcf_with_index)

    // -------------------------------------------------------------------------
    // Stage 5 – GENERATE_PHENOFILE
    // Builds a PLINK-format phenotype file from the resolved samplesheet.
    // Runs in parallel — no dependencies on downstream QC stages.
    // -------------------------------------------------------------------------
    GENERATE_PHENOFILE(PREP_INPUTS.out.resolved_samplesheet)

    // -------------------------------------------------------------------------
    // Stage 6 – MERGE_PLINK
    // Collects all per-plate bed/bim/fam files and merges them into a single
    // cohort-level PLINK dataset.
    // -------------------------------------------------------------------------
    ch_plink_beds = VCF_TO_PLINK.out.plink
        .map { plate, bed, bim, fam -> [bed, bim, fam] }
        .flatten()
        .collect()

    MERGE_PLINK(ch_plink_beds)

    // -------------------------------------------------------------------------
    // Stage 7 – SAMPLE_QC
    // Computes per-sample call rate and heterozygosity on the merged cohort.
    // Produces a keep-list of samples passing QC thresholds.
    // -------------------------------------------------------------------------
    SAMPLE_QC(MERGE_PLINK.out.merged)

    // -------------------------------------------------------------------------
    // Stage 8 – FILTER_SAMPLES
    // Removes samples failing SAMPLE_QC from the merged cohort.
    // Output: sample-filtered cohort with ALL SNPs intact (no SNP QC yet).
    // This is the last shared input before the pipeline forks into two parallel
    // branches: sex inference (Stages 9-10) and SNP QC (Stage 11).
    // -------------------------------------------------------------------------
    FILTER_SAMPLES(MERGE_PLINK.out.merged, SAMPLE_QC.out.keep_list)

    // -------------------------------------------------------------------------
    // Stage 9 – CHECK_SEX (GTC computed_gender)
    // Compares Illumina gender inference from GTC assessment with collected sex.
    // Runs after GTC_QC; does not require PLINK chrX extraction.
    // -------------------------------------------------------------------------
    CHECK_SEX(
        FILTER_GTC_SAMPLES.out.summary,
        ch_sex_info,
        PREP_INPUTS.out.resolved_samplesheet
    )

    // =========================================================================
    // ── SNP QC (autosomes only) ──────────────────────────────────────────────
    // =========================================================================

    // -------------------------------------------------------------------------
    // Stage 11 – SNP_QC
    // Applies --geno / --maf / --hwe to autosomes only.
    // Also removes duplicate SNP positions (--list-duplicate-vars).
    // -------------------------------------------------------------------------
    SNP_QC(FILTER_SAMPLES.out.plink)

    // =========================================================================
    // ── Downstream QC ────────────────────────────────────────────────────────
    // =========================================================================

    // -------------------------------------------------------------------------
    // Stage 12 – GENOTYPE_CONCORDANCE
    // All-vs-all pairwise genotype concordance on the SNP-QC-passed autosomal
    // cohort.  Uses chunked .traw streaming (bin/pairwise_concordance.py) so
    // peak memory stays flat regardless of cohort size.
    // Pairs >= params.concordance_dup_thresh flagged as likely duplicates.
    // -------------------------------------------------------------------------
    GENOTYPE_CONCORDANCE(SNP_QC.out.plink, file(params.samplesheet))

    // -------------------------------------------------------------------------
    // Stage 13 – IBD
    // Runs PLINK --genome on the merged cohort (pre-SNP-QC) after applying the
    // SAMPLE_QC keep-list.  Uses the merged dataset rather than the SNP-QC-
    // passed dataset so that IBD estimation includes all variants.
    // -------------------------------------------------------------------------
    IBD(
        MERGE_PLINK.out.merged,
        SAMPLE_QC.out.keep_list
    )

    // -------------------------------------------------------------------------
    // Stage 14 – PCA
    // Runs PLINK2 --pca on the merged cohort with the SAMPLE_QC keep-list.
    // Same reasoning as IBD: uses pre-SNP-QC merged data for full variant set.
    // -------------------------------------------------------------------------
    PCA(
        MERGE_PLINK.out.merged,
        SAMPLE_QC.out.keep_list
    )

    // -------------------------------------------------------------------------
    // Stage 15 – REPORT
    // Compiles all QC outputs into an HTML report + flagged-sample CSV.
    //
    // Inputs:
    //   sexcheck             – GTC computed_gender vs collected sex
    //   multimind            – per-sample discordance class
    //   plate_report         – per-plate discordance summary
    //   xy_tsv               – X/Y median intensities from GTC stats
    //   qc_stats             – per-sample call rate + het from SAMPLE_QC
    //   genome               – IBD pairwise estimates
    //   eigenvec             – PCA eigenvectors
    //   sex_info             – collected sex from samplesheet
    //   concordance          – pairwise genotype concordance results
    //   samplesheet          – original samplesheet for sample metadata lookup
    // -------------------------------------------------------------------------
    REPORT(
        CHECK_SEX.out.sexcheck,
        CHECK_SEX.out.multimind,
        CHECK_SEX.out.plate_report,
        CHECK_SEX.out.xy_tsv,
        SAMPLE_QC.out.qc_stats,
        IBD.out.genome,
        PCA.out.eigenvec,
        ch_sex_info,
        GENOTYPE_CONCORDANCE.out.concordance,
        file(params.samplesheet),
        FILTER_GTC_SAMPLES.out.summary,
        FILTER_GTC_SAMPLES.out.poor_gc10_list
    )
}