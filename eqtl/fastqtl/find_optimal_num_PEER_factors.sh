#!/bin/bash
# find optimal number of PEER factors

# variables: 
fastqtl=/usr/bin/fastQTL

# make appropriate directories: 
mkdir -p ../processed_data/eqtl/fastqtl/find_optimum_num_PEER_factors/


n=0
for num_geno_pc in 3 4 5; do
	for num_peer_factor in {1..15}; do
		# set covariate file name:
		covariates=../processed_data/160805/find_optimum_num_PEER_factors/covariates.fastqtl.pc$num_geno_pc.peer$num_peer_factor.tsv
		
		# create covariate file: 
		Rscript ../scripts/eqtl/utils/combine_covariates.R \
			--genotype_pc=../processed_data/genotype/genotype_pc/genotype_pcs.52samples.tsv \
			--peer=../processed_data/eqtl/peer/factors.tsv \
			--sample_info=../data/sample_info/sample_info.xlsx \
			--output=$covariates \
			--gender_coding=letter \
			--num_geno_pc=$num_geno_pc \
			--num_peer_factor=$num_peer_factor
		bgzip $covariates


		# inputs to fastqtl: 
		vcf=../data/joint3/recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf.id.gz
		bed=../processed_data/eqtl/fastqtl/expression/combined.filter.norm.bed.gz
		out=../processed_data/eqtl/fastqtl/find_optimum_num_PEER_factors/hcasmc.eqtl.pc$num_geno_pc.peer$num_peer_factor.chr20.txt
		cov=$covariates.gz
		region=chr20
		log=../processed_data/eqtl/fastqtl/find_optimum_num_PEER_factors/hcasmc.eqtl.pc$num_geno_pc.peer$num_peer_factor.chr20.log


		# run fastqtl: 
		$fastqtl --vcf $vcf \
			--bed $bed \
			--out $out \
			--cov $cov \
			--region $region > $log &



		# pause after launching 10 jobs:  
		n=$(($n+1))
		if [[ $n -ge 11 ]]; then
			wait 
			n=0
		fi


	done
done


n=0
for num_geno_pc in 3 4 5; do
	for num_peer_factor in {1..15}; do
		
		out=../processed_data/eqtl/fastqtl/find_optimum_num_PEER_factors/hcasmc.eqtl.pc$num_geno_pc.peer$num_peer_factor.chr20.txt

		# p-value adjustment: 
		Rscript ../scripts/eqtl/fastqtl/fastqtl_pvalue_corrections.nominal_pass.R $out ${out/txt/padj.txt} &

		# pause after launching 10 jobs:  
		n=$(($n+1))
		if [[ $n -ge 11 ]]; then
			wait 
			n=0
		fi
	done 
done 
