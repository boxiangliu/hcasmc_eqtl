#!/bin/bash
# find optimal number of PEER factors

# variables: 
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts
processed_data=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/
data=/srv/persistent/bliu2/HCASMC_eQTL/data
fastqtl=/usr/bin/fastQTL

# make appropriate directories: 
if [[ ! -d $processed_data/160805/find_optimum_num_PEER_factors ]]; then mkdir $processed_data/160805/find_optimum_num_PEER_factors; fi


n=0
for num_geno_pc in 3 4 5; do
	for num_peer_factor in {1..15}; do
		# set covariate file name:
		covariates=$processed_data/160805/find_optimum_num_PEER_factors/covariates.fastqtl.pc$num_geno_pc.peer$num_peer_factor.tsv
		
		# create covariate file: 
		Rscript $scripts/160530/combine_covariates.R \
			--genotype_pc=$processed_data/160519_genotype_PCA/genotype_pcs.52samples.tsv \
			--peer=$processed_data/160527/factors.tsv \
			--sample_info=$data/sample_info/sample_info.xlsx \
			--output=$covariates \
			--gender_coding=letter \
			--num_geno_pc=$num_geno_pc \
			--num_peer_factor=$num_peer_factor
		bgzip $covariates


		# inputs to fastqtl: 
		vcf=$data/joint3/recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf.id.gz
		bed=$processed_data/160530/combined.filter.norm.bed.gz
		out=$processed_data/160805/find_optimum_num_PEER_factors/hcasmc.eqtl.pc$num_geno_pc.peer$num_peer_factor.chr20.txt
		cov=$covariates.gz
		region=chr20
		log=$processed_data/160805/find_optimum_num_PEER_factors/hcasmc.eqtl.pc$num_geno_pc.peer$num_peer_factor.chr20.log
		

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
		
		out=$processed_data/160805/find_optimum_num_PEER_factors/hcasmc.eqtl.pc$num_geno_pc.peer$num_peer_factor.chr20.txt

		# p-value adjustment: 
		Rscript $scripts/160629/fastqtl_nominal_pvalue_corrections.R $out ${out/txt/padj.txt} &

		# pause after launching 10 jobs:  
		n=$(($n+1))
		if [[ $n -ge 11 ]]; then
			wait 
			n=0
		fi
	done 
done 
