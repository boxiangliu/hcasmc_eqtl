#!/bin/bash 
# bosh liu
# 2016/04/10

bed=$1
genetic_map=$2
reference_hap=$3
reference_legend=$4
reference_sample=$5
phased=$6
exclude=$7
shapeit -B ${bed/.bed/} -M $genetic_map -R $reference_hap $reference_legend $reference_sample -O ${phased/.haps/} --exclude-snp $exclude --thread 12 --window 0.5 --output-log $phased.log
