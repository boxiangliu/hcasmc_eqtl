for i in {1..22}; do 
	n_lines=$(wc -l ../processed_data/rasqual/input/rasqual.input.chr$i.txt | cut -d" " -f1)
	echo INFO - chr$i
	parallel -j10 bash rasqual/rasqual.sh \
	../processed_data/rasqual/input/rasqual.input.chr$i.txt {} \
	../processed_data/rasqual/expression/Y.tidy.bin \
	../processed_data/rasqual/expression/K.bin \
	../processed_data/rasqual/expression/X.bin \
	../data/joint3/asvcf/phased_and_imputed.chr$i.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz \
	../processed_data/rasqual/output/chr$i/ ::: `seq $n_lines`
done