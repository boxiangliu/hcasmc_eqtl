#!/bin/bash 
# bosh liu
# durga
# run beagle comform-gt 

# input arguments:
ref=$1
gt=$2
chrom=$3
out_prefix=$4
excludesamples=$5
match=POS
strict=false


# paths: 
conform_gt=/srv/persistent/bliu2/tools/beagle/conform-gt.24May16.cee.jar
java=/srv/persistent/bliu2/tools/jre1.8.0_91/bin/java


# run: 
$java -Xmx4g -jar $conform_gt \
	ref=$ref \
	gt=$gt \
	chrom=$chrom \
	out=$out_prefix \
	excludesamples=$excludesamples \
	strict=$strict \
	match=$match
