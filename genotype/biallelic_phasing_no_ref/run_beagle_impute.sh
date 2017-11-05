#!/bin/bash

input=$1
output_prefix=$2
# input=$wd/recalibrated_biallelic_variants.pass.vcf
# output_prefix=$wd/recalibrated_biallelic_variants.beagle

java=/srv/persistent/bliu2/tools/jre1.8.0_91/bin/java
beagle=/srv/persistent/bliu2/tools/beagle/beagle.03May16.862.jar

n=0
for i in $(seq 13 22);do
	n=$((n+1))
	if [[ n -gt 10 ]];then
		wait
		n=0
	fi 
	$java -Xmx4096m -jar $beagle nthreads=2 chrom=chr$i gl=$input out=$output_prefix.chr$i &
done
wait

# index each vcf:
for i in $(seq 1 22); do
tabix -p vcf $output_prefix.chr$i.vcf.gz
done 

# concatenate all chromosomes:
bcftools concat -Ov -o $output_prefix.vcf $output_prefix.chr{1..22}.vcf.gz

# remove intermediate files to folder: 
rm $output_prefix.chr{1..22}.vcf.gz
rm $output_prefix.chr{1..22}.vcf.gz.tbi
rm $output_prefix.chr{1..22}.log
