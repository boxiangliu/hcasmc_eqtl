# prepare covariates: 
# top three genotype PCs
# 15 PEER factors 
# gender
# example: /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_data/eQTLInputFiles/covariates/
Rscript eqtl/utils/combine_covariates.R \
	--genotype_pc=../processed_data/genotype/genotype_pc/genotype_pcs.52samples.tsv \
	--peer=../processed_data/eqtl/peer/factors.tsv \
	--sample_info=/srv/persistent/bliu2/HCASMC_eQTL/data/sample_info/sample_info.xlsx \
	--output=../processed_data/eqtl/matrix_eqtl/combine_covariates/covariates.tsv \
	--gender_coding=numerical