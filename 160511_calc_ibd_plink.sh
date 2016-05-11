#!/bin/bash 
# 2016/05/11
# bosh liu
# durga
# calculate IBD using plink1.9

# paths: 
vcf=../data/joint/recalibrated_variants.GRCh37.pass.vcf.gz
output_dir=../processed_data/160511_calc_ibd_plink

# select independent SNPs:
if [[ ! -e $output_dir/indep_pairwise_50_5_0.2.prune.in ]]; then
	plink --vcf $vcf --keep-allele-order --indep-pairwise 50 5 0.2 --out $output_dir/indep_pairwise_50_5_0.2
fi 

# create independent SNP set: 
if [[ ! -e $output_dir/indep_pairwise_50_5_0.2.bim ]]; then 
	plink --vcf $vcf --make-bed --extract $output_dir/indep_pairwise_50_5_0.2.prune.in --out $output_dir/indep_pairwise_50_5_0.2
fi

# calculate IBS and IBD:
plink --bfile $output_dir/indep_pairwise_50_5_0.2 --genome --out $output_dir/indep_pairwise_50_5_0.2