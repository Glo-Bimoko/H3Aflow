# Archived pipeline code

These files supported the **PLINK chrX sex-check branch** (`SPLIT_CHROM` → stratified chrX QC → `--check-sex` → multimind), replaced by **`bin/gtc_sex_check.py`** using Illumina GTC `computed_gender` from `GTC_QC` stats.

## `bin/`

| Script | Role |
|--------|------|
| `filter_chrx_snps.py` | Sex-stratified chrX SNP QC before `--check-sex` |
| `make_sex_update.py` | Built PLINK `--update-sex` file from `sex_info.tsv` |
| `rename_chr23.py` | Renamed chr 23 → X in BIM for `--split-x` |
| `sexcheck_multimind.py` | Multi-`--mind` `--check-sex` (H3AGWAS xCheck style) |
| `sexcheck_plate_report.py` | Per-plate discordance from PLINK sexcheck |
| `annotate_sex_check.py` | Merged PLINK sexcheck with collected sex |
| `compare_sex.py` | F-statistic plots vs collected sex |
| `extract_xy_intensity.py` | X/Y intensities from concatenated GTC TSVs |
| `write_sexcheck_stubs.py` | Empty outputs when chrX check could not run |
| `write_placeholder_plot.py` | Placeholder PNG when plotting failed |

## `modules/`

| Module | Role |
|--------|------|
| `split_chrom.nf` | Extracted chrX/chrY PLINK subsets |
| `xy_intensity.nf` | Ran `extract_xy_intensity.py` on GTC TSVs |

Active sex check: `modules/check_sex.nf` → `bin/gtc_sex_check.py`.
