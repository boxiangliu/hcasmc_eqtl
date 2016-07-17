#!/bin/bash 
htseq=/srv/persistent/bliu2/tools/HTSeq-0.6.1/scripts
# annotation=/srv/persistent/bliu2/shared/annotation/gencode/v14/gencode.v14.annotation.gtf
annotation=/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf

input_dir=$1
input_files=$2
output_dir_prefix=$3
# input_dir=../data/rnaseq/alignments/
# output_dir_prefix=../data/rnaseq/expression/
num_processes=0

while read input; do
	if [[ $input == \#* ]]; then continue; fi
	echo "input:" $input_dir/$input

	# forward:
	output_dir=$output_dir_prefix/forward
	[[ ! -d $output_dir ]] && mkdir -p $output_dir

	# gene level quantification: 
	output=${input/bam/gene.count}
	echo "output:" $output_dir/$output
	[[ ! -d $(dirname $output_dir/$output) ]] && mkdir -p $(dirname $output_dir/$output)
	($htseq/htseq-count --format=bam --order=pos --mode=intersection-strict --stranded=yes --idattr=gene_id --type=gene --order=pos $input_dir/$input $annotation > $output_dir/$output) &
	num_processes=$((num_processes+1))

	# exon level quantification: 
	# output=${input/bam/exon.count}
	# echo "output:" $output_dir/$output
	# [[ ! -d $(dirname $output_dir/$output) ]] && mkdir -p $(dirname $output_dir/$output)
	# ($htseq/htseq-count --format=bam --order=pos --mode=intersection-strict --stranded=yes --idattr=gene_id --type=exon --order=pos $input_dir/$input $annotation > $output_dir/$output) &
	# num_processes=$((num_processes+1))

	# reverse:
	output_dir=$output_dir_prefix/reverse
	[[ ! -d $output_dir ]] && mkdir -p $output_dir

	# gene level quantification: 
	output=${input/bam/gene.count}
	echo "output:" $output_dir/$output
	[[ ! -d $(dirname $output_dir/$output) ]] && mkdir -p $(dirname $output_dir/$output)
	($htseq/htseq-count --format=bam --order=pos --mode=intersection-strict --stranded=reverse --idattr=gene_id --type=gene --order=pos $input_dir/$input $annotation > $output_dir/$output) &
	num_processes=$((num_processes+1))


	# exon level quantification: 
	# output=${input/bam/exon.count}
	# echo "output:" $output_dir/$output
	# [[ ! -d $(dirname $output_dir/$output) ]] && mkdir -p $(dirname $output_dir/$output)
	# ($htseq/htseq-count --format=bam --order=pos --mode=intersection-strict --stranded=reverse --idattr=gene_id --type=exon --order=pos $input_dir/$input $annotation > $output_dir/$output) &
	# num_processes=$((num_processes+1))


	if [[ $num_processes -ge 10 ]]; then
		wait 
		num_processes=0
	fi 

done < $input_files
