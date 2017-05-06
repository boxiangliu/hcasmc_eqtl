#!/bin/bash 
# bosh liu
# durga
# run beagle phasing without reference

input=$1
out_prefix=$2
map=$3
chrom=$4
impute=false
beagle=/srv/persistent/bliu2/tools/beagle/beagle.03May16.862.jar
java=/srv/persistent/bliu2/tools/jre1.8.0_91/bin/java


$java -Xmx4096m -jar $beagle \
	nthreads=1 \
	chrom=$chrom \
	gt=$input \
	out=$out_prefix \
	map=$map \
	impute=$impute

