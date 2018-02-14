# eQTL analysis pipeline
# 2017-11-22
# Boxiang Liu

# Setup:
mkdir eqtl ../processed_data/eqtl ../figures/eqtl

# extract PEER factors:
bash eqtl/peer/peer_factor.sh
Rscript eqtl/peer/peer_factor_correlation.R


# Matrix eQTL:
bash eqtl/matrix_eqtl/combine_covariates.sh
Rscript eqtl/matrix_eqtl/covariates_correlation.R
bash eqtl/matrix_eqtl/matrix_eqtl.sh


# find the best combination of genotype PCs and PEER factors:
bash eqtl/fastqtl/find_optimal_num_PEER_factors.sh 
Rscript eqtl/fastqtl/plot_num_egene_vs_cov.R


# create covariate file: 
Rscript eqtl/utils/combine_covariates.R \
	--genotype_pc=../processed_data/genotype/genotype_pc/genotype_pcs.52samples.tsv \
	--peer=../processed_data/eqtl/peer/factors.tsv \
	--sample_info=../data/sample_info/sample_info.xlsx \
	--output=../processed_data/eqtl/fastqtl/covariates/covariates.tsv \
	--gender_coding=letter \
	--num_geno_pc=4 \
	--num_peer_factor=8
bgzip ../processed_data/eqtl/fastqtl/covariates/covariates.tsv


# FastQTL:
bash eqtl/fastqtl/change_sid.sh
bash eqtl/fastqtl/fastqtl.nominal.sh
bash eqtl/fastqtl/adjust_pvalue.nominal.sh

# FastQTL permutation:
bash eqtl/fastqtl/fastqtl.perm.sh 
Rscript eqtl/fastqtl/count_sig_eqtl.sh

# plot the p-values: 
Rscript sqtl/fastQTL/plot_pvalue.R \
../processed_data/eqtl/fastqtl/output/nominal/all.txt.gz \
../figures/eqtl/fastqtl/plot_pvalue/

# plot the number of eQTLs vs distance to TSS: 
Rscript eqtl/fastqtl/plot_eqtl_vs_distance.R \
	-eqtl_file=../processed_data/eqtl/fastqtl/output/nominal/all.txt.gz \
	-pval_vs_dist=../figures/eqtl/fastqtl/plot_eqtl_vs_distance/eqtl_pval_vs_dist.pdf

# Enrichment analysis:
bash eqtl/enrichment/vsea.R