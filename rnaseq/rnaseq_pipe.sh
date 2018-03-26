# Setup: 
mkdir ../data/rnaseq2/
mkdir ../data/rnaseq_dase/

# Preprocessing (for eQTL samples):
bash rnaseq/preprocess/get_rnaseq_data.sh
bash rnaseq/preprocess/sort_and_index.sh
bash rnaseq/preprocess/quickcheck.sh
bash rnaseq/preprocess/calculate_median_coverage.sh
bash rnaseq/preprocess/RNAseQC.sh
bash rnaseq/preprocess/copy_rpkm.sh
bash rnaseq/preprocess/combine_rpkm.sh
Rscript rnaseq/preprocess/filter_rpkm.R 0.1 10 ../processed_data/rnaseq/preprocess/combine_rpkm/combined.rpkm ../processed_data/rnaseq/preprocess/combine_rpkm/combined.filter.rpkm


# Preprocessing (for dASE samples):
bash rnaseq/preprocess_dASE/copy_rnaseq_from_valk.sh # copy dASE RNAseq from valk
bash rnaseq/preprocess_dASE/rename_rnaseq_samples.sh # rename dASE samples
bash rnaseq/preprocess_dASE/RNAseQC.sh # run RNAseQC
bash rnaseq/preprocess_dASE/copy_rpkm.sh # copy rpkm to the rpkm folder
bash rnaseq/preprocess_dASE/combine_rpkm.sh # combine FBS and SF rpkm to two separate files


# Quality control:
Rscript rnaseq/quality_control/sample_correlation.R
