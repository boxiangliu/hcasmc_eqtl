#!/bin/bash 
# bosh liu
# durga
# 2016/05/11
# calculate allele frequency of the pruned marker set


# paths: 
vcf=../data/joint/recalibrated_variants.GRCh37.pass.vcf.gz
input_dir=../processed_data/160511_calc_ibd_plink
output_dir=../processed_data/160511_calc_freq


# calculate frequency:
plink --vcf $vcf --freq --extract $input_dir/indep_pairwise_50_5_0.2.prune.in --out $output_dir/indep_pairwise_50_5_0.2