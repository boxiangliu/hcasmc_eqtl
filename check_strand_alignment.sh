#!/bin/bash
# bosh liu
# 2016/04/10

bed=$1
reference_hap=$2
reference_legend=$3
reference_sample=$4
output=$5

shapeit -check --input-bed ${bed/.bed/} --input-ref $reference_hap $reference_legend $reference_sample --output-log ${output/.snp.strand/} > $bed.check_strand_alignment.log
