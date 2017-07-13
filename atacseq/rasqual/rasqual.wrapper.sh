in_dir=../processed_data/atacseq/rasqual/

for i in {1..22}; do 
	echo INFO - chr$i
	mkdir -p $in_dir/output/chr$i/

	n_lines=$(wc -l $in_dir/input/rasqual.input.chr$i.txt | cut -d" " -f1)

	parallel -j15 bash atacseq/rasqual/rasqual.sh \
	$in_dir/input/rasqual.input.chr$i.txt {} \
	$in_dir/expression/Y.bin \
	$in_dir/expression/K.bin \
	../data/joint3/asvcf_sid_atac/phased_and_imputed.chr$i.rename.dr2.hwe.indellt51.atacsample.hg19.vcf.new.gz \
	$in_dir/output/chr$i/ ::: `seq $n_lines`
done