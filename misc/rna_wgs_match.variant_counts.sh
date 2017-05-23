#!/bin/bash 
# bosh liu
# 2016/04/06
# convert mpile to counts 

# it is normal to see the following errors:
# Traceback (most recent call last):
#   File "mpileupCountsTaregtedSNPs.py", line 63, in <module>
#     line=inFile.next()
# StopIteration
# Or
# Traceback (most recent call last):
#   File "mpileupCountsTaregtedSNPs.py", line 13, in <module>
#     line=inFile.next()
# StopIteration

MPILEUP_DIR=$1
SAMPLE_LIST=$2
VARIANT_COUNT_DIR=$3
while read sample; do
	[[ ! -d $(dirname $VARIANT_COUNT_DIR/$sample) ]] && mkdir -p $(dirname $VARIANT_COUNT_DIR/$sample)
	echo $sample
	python mpileupCountsTaregtedSNPs.py $MPILEUP_DIR/$sample > $VARIANT_COUNT_DIR/${sample/pileup/count}
done < $SAMPLE_LIST
touch .rna_wgs_match.variant_counts.sh.done