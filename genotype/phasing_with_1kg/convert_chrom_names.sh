#!/bin/bash
# bosh liu
# durga
# convert chromsome names
# e.g. from hg19 to GRCh37

# command args: 
input=$1
output=$2
old_to_new=$3

bcftools annotate --rename-chrs $old_to_new -Oz -o $output $input
