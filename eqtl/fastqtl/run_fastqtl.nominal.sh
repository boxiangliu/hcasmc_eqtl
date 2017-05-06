#!/usr/bin/env Rscript 
# bosh liu
# durga
# run fastqtl

vcf=$1
bed=$2
cov=$3
out=$4
region=$5

# fastqtl=/srv/persistent/bliu2/tools/FastQTL/bin/fastQTL.static
fastqtl=/usr/bin/fastQTL
processed_data=../processed_data/160530/


# $fastqtl --vcf /srv/persistent/bliu2/HCASMC_eQTL/data/joint2/recalibrated_biallelic_SNP.beagle.rename.dr2.hwe.maf.vcf.gz \
# 	--bed $processed_data/combined.filter.norm.bed.gz \
# 	--out $processed_data/permutation.txt.gz \
# 	--cov $processed_data/covariates.fastqtl.tsv.gz \
# 	--permute 1000 10000 \
# 	--region chr22

echo "$fastqtl --vcf $vcf \
	--bed $bed \
	--out $out \
	--cov $cov \
	--region $region"

$fastqtl --vcf $vcf \
	--bed $bed \
	--out $out \
	--cov $cov \
	--normal \
	--region $region

