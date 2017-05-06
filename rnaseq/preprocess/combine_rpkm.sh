#!/bin/bash 
src=/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/rpkm
dst=../processed_data/rnaseq/preprocess/combine_rpkm/

tail -n +3 $src/1020301/genes.rpkm | cut -f1 > $dst/combined.rpkm
samples=($(ls -d $src/*/))
for sample in ${samples[@]};do
	sample=$(basename $sample)
	sample=${sample///}
	echo $sample
	tail -n +3 $src/$sample/genes.rpkm | cut -f3 > $dst/$sample.rpkm.tmp
	cp $dst/combined.rpkm $dst/combined.rpkm.tmp
	paste -d "\t" $dst/combined.rpkm.tmp $dst/$sample.rpkm.tmp > $dst/combined.rpkm
done
rm $dst/*.tmp