#!/bin/bash 
# bosh liu
# durga
# beagle phasing with reference

input=$1
out_prefix=$2
ref=$3
map=$4
chrom=$5
impute=false
beagle=/srv/persistent/bliu2/tools/beagle/beagle.03May16.862.jar
java=/srv/persistent/bliu2/tools/jre1.8.0_91/bin/java


$java -Xmx4096m -jar $beagle \
	nthreads=2 \
	chrom=$chrom \
	gt=$input \
	out=$out_prefix \
	ref=$ref \
	map=$map \
	impute=$impute

