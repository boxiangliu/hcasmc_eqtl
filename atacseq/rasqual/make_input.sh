Rscript atacseq/rasqual/gtf.R


out_dir=../processed_data/atacseq/rasqual/
mkdir -p $out_dir/input $out_dir/expression/

# Prepare input: 
parallel -j10 grep chr{} \
../processed_data/atacseq/rasqual/merged_peak.gtf "|" \
python rasqual/make_input.py \
../data/joint3/asvcf_sid_atac/phased_and_imputed.chr{}.rename.dr2.hwe.indellt51.atacsample.hg19.vcf.new.gz 100000 '>' \
$out_dir/input/rasqual.input.chr{}.txt ::: {1..22} X


cat $out_dir/input/rasqual.input.chr{1..22}.txt > $out_dir/input/rasqual.input.autosome.txt


# Calculate GC concent: 
python rasqual/calc_gcc.py \
/mnt/lab_data/montgomery/shared/genomes/hg19/hg19.fa \
../processed_data/atacseq/rasqual/merged_peak.gtf \
exon > $out_dir/expression/gcc.txt


# Make count matrix: 
cat ../data/atacseq/fbs/count/1346.count | awk '{print $1"_"$2"_"$3}' > $out_dir/colname.txt
parallel cut -f4 ../data/atacseq/fbs/count/{}.count '>' $out_dir/{}.txt ::: 1346  1508  1522  200212  2305  2356  2510  2989
cd $out_dir
paste colname.txt 1346.txt  1508.txt  1522.txt  200212.txt  2305.txt  2356.txt  2510.txt  2989.txt > expression/Y.txt
rm {1346,1508,1522,200212,2305,2356,2510,2989}.txt
rm colname.txt


# Calculate RASQUAL offset:
cut -f2-2 $out_dir/expression/gcc.txt > $out_dir/expression/gcc.cut.txt
R --vanilla --quiet --args $out_dir/expression/Y.txt $out_dir/expression/gcc.cut.txt $out_dir/expression/K.txt < rasqual/makeOffset.R


# Binarize: 
R --vanilla --quiet --args $out_dir/expression/Y.txt $out_dir/expression/K.txt < ../../tools/rasqual/R/txt2bin.R 
