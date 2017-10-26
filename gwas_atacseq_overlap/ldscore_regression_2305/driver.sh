# Using 2305 only insteading of merging all 9 ATACseq samples:
# bash gwas_atacseq_overlap/ldscore_regression_2305/merge_peaks.R


# Make tissue specific snp annotation: 
# Rscript gwas_atacseq_overlap/ldscore_regression/tissue_specific_snp_annotation.R ../processed_data/gwas_atacseq_overlap/ldscore_regression_2305/merge_peaks_released_all_cell_type/ ../processed_data/gwas_atacseq_overlap/ldscore_regression_2305/tissue_specific_snp_annotation/


# Estimating LD score:
bash gwas_atacseq_overlap/ldscore_regression/ldscore_merged.sh ../processed_data/gwas_atacseq_overlap/ldscore_regression_2305/tissue_specific_snp_annotation/ ../processed_data/gwas_atacseq_overlap/ldscore_regression_2305/ldscore_merged/ ../logs/gwas_atacseq_overlap/ldscore_regression_2305/ldscore_merged/


# Estimate heritability:
bash gwas_atacseq_overlap/ldscore_regression/partition_heritability_merged.sh ../processed_data/gwas_atacseq_overlap/ldscore_regression_2305/ldscore_merged/ ../processed_data/gwas_atacseq_overlap/ldscore_regression_2305/partition_heritability_merged/
