#!/bin/bash
# find optimal number of PEER factors

scripts=./160530
processed_data=../processed_data/160530
n=0
for num_geno_pc in 3 4 5; do
	for num_peer_factor in {1..15}; do
		covariates=$processed_data/find_optimum_num_PEER_factors/covariates.fastqtl.pc$num_geno_pc.peer$num_peer_factor.tsv
		Rscript $scripts/combine_covariates.R \
			--genotype_pc=../processed_data/160519_genotype_PCA/genotype_pcs.52samples.tsv \
			--peer=../processed_data/160527/factors.tsv \
			--sample_info=/srv/persistent/bliu2/HCASMC_eQTL/data/sample_info/sample_info.xlsx \
			--output=$covariates \
			--gender_coding=letter \
			--num_geno_pc=$num_geno_pc \
			--num_peer_factor=$num_peer_factor
		bgzip $covariates

		for i in {1..22}; do 
			bash $scripts/run_fastqtl.sh \
				/srv/persistent/bliu2/HCASMC_eQTL/data/joint2/recalibrated_biallelic_SNP.beagle.rename.dr2.hwe.maf.vcf.id.gz \
				$processed_data/combined.filter.norm.bed.gz \
				$covariates.gz \
				$processed_data/find_optimum_num_PEER_factors/fastqtl.pc$num_geno_pc.peer$num_peer_factor.chr$i.txt.gz \
				chr$i \
				1000 \
				10000 > $processed_data/find_optimum_num_PEER_factors/run_fastqtl.pc$num_geno_pc.peer$num_peer_factor.chr$i.log &
			
			n=$(($n+1))
			if [[ $n -ge 11 ]]; then
				wait 
				n=0
			fi
		done
	done
done
wait

for num_geno_pc in 3 4 5; do
	for num_peer_factor in {1..15}; do
		zcat $processed_data/find_optimum_num_PEER_factors/fastqtl.pc$num_geno_pc.peer$num_peer_factor.chr{1..22}.txt.gz > $processed_data/find_optimum_num_PEER_factors/fastqtl.pc$num_geno_pc.peer$num_peer_factor.txt
		Rscript $scripts/fastqtl_pvalue_corrections.R $processed_data/find_optimum_num_PEER_factors/fastqtl.pc$num_geno_pc.peer$num_peer_factor.txt $processed_data/find_optimum_num_PEER_factors/fastqtl.pc$num_geno_pc.peer$num_peer_factor.padj.txt
	done 
done 