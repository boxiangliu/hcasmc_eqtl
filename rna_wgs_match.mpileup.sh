#!/bin/bash 
# Bosh Liu
# 2016/04/06
# run samtools mpileup on variant sites 

BAM_DIR=$1
SAMPLE_LIST=$2
PILEUP_DIR=$3
VARIANT=$4

QUAL=30 
REF=/srv/persistent/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

# create pileup directory if needed: 
[[ ! -d $PILEUP_DIR ]] && mkdir -p $PILEUP_DIR

while read sample; do
	echo $sample
	[[ ! -d  $(dirname $PILEUP_DIR/$sample.pileup) ]] && mkdir -p $(dirname $PILEUP_DIR/$sample.pileup)
	echo "writing to $PILEUP_DIR/$sample.pileup"
	samtools mpileup -Q $QUAL -BAf $REF -l $VARIANT $BAM_DIR/$sample > $PILEUP_DIR/${sample/out.bam/pileup}
done < $SAMPLE_LIST

touch .rna_wgs_match.mpileup.sh.done