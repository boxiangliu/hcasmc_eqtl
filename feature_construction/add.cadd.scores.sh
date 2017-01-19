#!/bin/bash

# Add CADD scores to each variant BED with MAFs and allele info

# takes as input the chromosome number
# outputs a bed file with CADD score (raw and phred) for each variant in the bed file for each chromosome
# if no CADD score is available, outputs NA


in_file=$1
out_file=$2
cat $in_file | \
python feature_construction/extractCADDscores.py - > $out_file
gzip $out_file



