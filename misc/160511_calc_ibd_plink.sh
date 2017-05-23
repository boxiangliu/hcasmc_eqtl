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
	plink --vcf $vcf --keep-allele-order --maf 0.05 --geno 0.1 --indep-pairwise 50 5 0.2 --out $output_dir/indep_pairwise_50_5_0.2
fi 

# remove "." in the prune.in file: 
awk '{ if ($1 != ".") print $0}' $output_dir/indep_pairwise_50_5_0.2.prune.in > $output_dir/indep_pairwise_50_5_0.2.prune.in.nodot


# calculate IBS and IBD:
# plink --bfile $output_dir/indep_pairwise_50_5_0.2 --genome --out $output_dir/indep_pairwise_50_5_0.2
plink --vcf $vcf --extract $output_dir/indep_pairwise_50_5_0.2.prune.in.nodot --genome --out $output_dir/indep_pairwise_50_5_0.2
