#!/bin/bash
# find optimal number of PEER factors for matrix eQTL

scripts=./160530
processed_data=../processed_data/160530
n=0
# for num_geno_pc in 3 4 5; do
# 	for num_peer_factor in {1..15}; do
for num_geno_pc in 3 4 5; do
	for num_peer_factor in {1..15}; do
		covariates=$processed_data/find_optimum_num_PEER_factors_matrixeqtl/covariates.matrixeqtl.pc$num_geno_pc.peer$num_peer_factor.tsv
		# Rscript $scripts/combine_covariates.R \
		# 	--genotype_pc=../processed_data/160519_genotype_PCA/genotype_pcs.52samples.tsv \
		# 	--peer=../processed_data/160527/factors.tsv \
		# 	--sample_info=/srv/persistent/bliu2/HCASMC_eQTL/data/sample_info/sample_info.xlsx \
		# 	--output=$covariates \
		# 	--gender_coding=numerical \
		# 	--num_geno_pc=$num_geno_pc \
		# 	--num_peer_factor=$num_peer_factor

		# run matrix eqtl:
		# if [[ ! -s $processed_data/find_optimum_num_PEER_factors_matrixeqtl/pc$num_geno_pc.peer$num_peer_factor.cis.txt ]]; then
		Rscript $scripts/run_matrix_eQTL.R \
			$processed_data/dosage.tsv \
			$processed_data/genotype_loc.txt \
			$processed_data/combined.filter.norm.rpkm \
			$processed_data/gene_loc.txt \
			$covariates \
			$processed_data/find_optimum_num_PEER_factors_matrixeqtl/pc$num_geno_pc.peer$num_peer_factor.2.
		# fi


		# n=$(($n+1))
		# if [[ $n -ge 5 ]]; then
		# 	wait 
		# 	n=0
		# fi

	done
done
wait

for num_geno_pc in 3 4 5; do
	for num_peer_factor in {1..15}; do
		# plot histogram and qqplot of q-value: 
		Rscript $scripts/qqplot_pvalue.R \
			$processed_data/find_optimum_num_PEER_factors_matrixeqtl/pc$num_geno_pc.peer$num_peer_factor.2.cis.txt \
			$processed_data/find_optimum_num_PEER_factors_matrixeqtl/pc$num_geno_pc.peer$num_peer_factor.2
	done
done 