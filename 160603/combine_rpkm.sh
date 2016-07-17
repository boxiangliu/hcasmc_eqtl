#!/bin/bash 
src=$1
dst=$2
src=/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq_dase/rpkm
dst=../processed_data/160603/rpkm

# fbs samples:
tail -n +3 $src/1051601_fbs/genes.rpkm | cut -f1 > $dst/dase.fbs.rpkm
samples=($(ls -d $src/*_fbs*/))
for sample in ${samples[@]};do
	sample=$(basename $sample)
	sample=${sample///}
	echo $sample
	tail -n +3 $src/$sample/genes.rpkm | cut -f3 > $dst/$sample.rpkm.tmp
	cp $dst/dase.fbs.rpkm $dst/dase.fbs.rpkm.tmp
	paste -d "\t" $dst/dase.fbs.rpkm.tmp $dst/$sample.rpkm.tmp > $dst/dase.fbs.rpkm
done
rm $dst/*.tmp


# sf samples: 
tail -n +3 $src/1051601_sf/genes.rpkm | cut -f1 > $dst/dase.sf.rpkm
samples=($(ls -d $src/*_sf*/))
for sample in ${samples[@]};do
	sample=$(basename $sample)
	sample=${sample///}
	echo $sample
	tail -n +3 $src/$sample/genes.rpkm | cut -f3 > $dst/$sample.rpkm.tmp
	cp $dst/dase.sf.rpkm $dst/dase.sf.rpkm.tmp
	paste -d "\t" $dst/dase.sf.rpkm.tmp $dst/$sample.rpkm.tmp > $dst/dase.sf.rpkm
done
rm $dst/*.tmp