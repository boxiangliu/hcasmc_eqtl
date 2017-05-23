#!/bin/bash 
# bosh liu
# 2016/04/10

input=$1
bed=$2
chr=$3

plink \
--vcf $input \
--keep-allele-order \
--make-bed \
--freq \
--chr $chr \
--out ${bed/.bed/}