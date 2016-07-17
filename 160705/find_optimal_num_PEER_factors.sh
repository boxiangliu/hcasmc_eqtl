#!/bin/bash
# find optimal number of PEER factors
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts/
processed_data=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/
figures=/srv/persistent/bliu2/HCASMC_eQTL/figures/
for i in {3,6,9,12,15}; do
	# # create covariate file:
	# covariates=$processed_data/160705/covariates.top10000.peer$i.fastqtl.tsv
	# Rscript $scripts/160530/combine_covariates.R \
	# 	--genotype_pc=../processed_data/160519_genotype_PCA/genotype_pcs.52samples.tsv \
	# 	--peer=$processed_data/160705/factors.tsv \
	# 	--sample_info=/srv/persistent/bliu2/HCASMC_eQTL/data/sample_info/sample_info.xlsx \
	# 	--output=$covariates \
	# 	--gender_coding=letter \
	# 	--num_geno_pc=3 \
	# 	--num_peer_factor=$i
	# bgzip $covariates


	# # run fastqtl:
	# bash $scripts/160629/run_fastqtl.nominal.wrap.sh \
	# /srv/persistent/bliu2/HCASMC_eQTL/data/joint3/recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf.id.gz \
	# $processed_data/160705/leafcutter.norm.perm.bed.gz \
	# $covariates.gz \
	# $processed_data/160705/sqtl.nominal.allpairs.normal.1e5.perm.top10000.peer$i.txt \
	# "--normal --window 1e5"


	# make qqplot and histograms:
	Rscript $scripts/160629/qqplot_pvalue.R \
	$processed_data/160705/sqtl.nominal.allpairs.normal.1e5.perm.top10000.peer$i.txt \
	$figures/160705/sqtl.perm.qqplot.top10000.peer$i.png \
	$figures/160705/sqtl.perm.histogram.top10000.peer$i.pdf
done