#!/bin/bash 
wd=$1
echo "working directory:" $wd
input=$2
echo "input:" $input
output=$3
echo "output:" $output
cd $wd 
HG19=/srv/persistent/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
CHROM_SIZE=/srv/persistent/bliu2/shared/genomes/hg19.chrom.sizes
while read sample; do
	sample=${sample}/recal_reads.bam
	bamToBed -i ${sample} > ${sample/bam/bed}
	sort -k1,1 ${sample/bam/bed} > ${sample/bam/sorted.bed}
	genomeCoverageBed -bga -i ${sample/bam/sorted.bed} -g ${HG19} > ${sample/bam/bedGraph}
	/srv/persistent/bliu2/tools/ucsc_tools/bedGraphToBigWig ${sample/bam/bedGraph} ${CHROM_SIZE} ${sample/bam/bw}
done < ${input}
touch ${output}
