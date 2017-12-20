#!/bin/bash 
src=/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq_dase/alignments/
dst=/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq_dase/rpkm/


samples=($(ls -d $src/*/))
echo "${#samples[@]} in total"
for sample in ${samples[@]}; do 
	sample=$(basename $sample)
	sample=${sample///}
	echo $sample
	mkdir $dst/$sample
	cp $src/$sample/report/genes.rpkm.gct $dst/$sample/genes.rpkm
done 
