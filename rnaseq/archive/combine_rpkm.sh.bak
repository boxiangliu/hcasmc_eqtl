# combine all RPKMs:
wd=../data/rnaseq2/alignments

tail -n +3 $wd/1020301/report/genes.rpkm.gct | cut -f1-2 > ../processed_data/160519_rpkm/combined.rpkm
samples=($(ls -d $wd/*/))

for sample in ${samples[@]};do
	sample=$(basename $sample)
	sample=${sample///}
	echo $sample
	tail -n +3 $wd/$sample/report/genes.rpkm.gct | cut -f3 > ../processed_data/160519_rpkm/$sample.rpkm.tmp
	cp ../processed_data/160519_rpkm/combined.rpkm ../processed_data/160519_rpkm/combined.rpkm.tmp
	paste -d "\t" ../processed_data/160519_rpkm/combined.rpkm.tmp ../processed_data/160519_rpkm/$sample.rpkm.tmp > ../processed_data/160519_rpkm/combined.rpkm
done
rm ../processed_data/160519_rpkm/*.tmp
