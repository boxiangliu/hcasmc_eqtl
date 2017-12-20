#!/bin/bash
input=rnaseq/preprocess_dASE/rename_rnaseq_samples.map.txt
dir=/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq_dase/alignments/
while read line; do 
	line=($line)
	old=${line[0]}
	new=${line[1]}
	echo $old 
	echo $new 
	mv $dir/$old $dir/$new
done < $input