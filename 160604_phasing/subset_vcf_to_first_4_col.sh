#!/bin/bash 
# durga
# bosh liu
# subset the chrom, pos, ref, alt field in vcf file

vcf=$1
output=$2
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\n' $vcf > $output

