region=chr6:133210276-135210276
tss=134210276
n_sample=52
n_rsnp=11130
n_fsnp=20
gene_id=ENSG00000118526.6
gene_name=TCF21
feat_id=$(grep $gene_id -n /srv/persistent/bliu2/HCASMC_eQTL/processed_data//rasqual/Y.tidy.txt | cut -d":" -f1,1)
exon_start_positions=134210276,134213075
exon_end_positions=134210985,134216691
window_size=2000000
vcf_file=../processed_data/160604_phasing/phased_and_imputed/phased_and_imputed.chr6.rename.dr2.indellt51.rnasamples.hg19.vcf.new.gz
tabix $vcf_file $region | \
	/srv/persistent/bliu2/tools/rasqual/src/rasqual \
	-y ../processed_data/rasqual/Y.tidy.bin \
	-k ../processed_data/rasqual/K.bin \
	-x ../processed_data/rasqual/X.bin \
	-n 52 -j $feat_id -l $n_rsnp -m $n_fsnp \
    -s $exon_start_positions -e $exon_end_positions \
    --imputation-quality 0.8 --imputation-quality-fsnp 0.8 \
    -f $gene_name -z --n_threads 1 \
    -c $tss -w $window_size -v --force > ../processed_data/rasqual/output/${gene_id}_${gene_name}_w${window_size}.txt


# conditional analysis on rs12190287: 
window_size=100000
tabix $vcf_file $region | \
	/srv/persistent/bliu2/tools/rasqual/src/rasqual \
	-y ../processed_data/rasqual/Y.tidy.bin \
	-k ../processed_data/rasqual/K.bin \
	-x ../processed_data/rasqual/X.bin \
	-n 52 -j $feat_id -l $n_rsnp -m $n_fsnp \
    -s $exon_start_positions -e $exon_end_positions \
    --imputation-quality 0.8 --imputation-quality-fsnp 0.8 \
    -f $gene_name -z --n_threads 1 \
    -c $tss -w $window_size -v -k2 rs12190287:0.701354 > ../processed_data/rasqual/output/${gene_id}_${gene_name}_w${window_size}_crs12190287.txt


# conditional analysis on rs2327429: 
tabix $vcf_file $region | \
	/srv/persistent/bliu2/tools/rasqual/src/rasqual \
	-y ../processed_data/rasqual/Y.tidy.bin \
	-k ../processed_data/rasqual/K.bin \
	-x ../processed_data/rasqual/X.bin \
	-n 52 -j $feat_id -l $n_rsnp -m $n_fsnp \
    -s $exon_start_positions -e $exon_end_positions \
    --imputation-quality 0.8 --imputation-quality-fsnp 0.8 \
    -f $gene_name -z --n_threads 1 \
    -c $tss -w $window_size -v -k2 rs2327429:0.701665 > ../processed_data/rasqual/output/${gene_id}_${gene_name}_w${window_size}_crs2327429.txt



tabix $vcf_file $region | \
	/srv/persistent/bliu2/tools/rasqual/src/rasqual \
	-y ../processed_data/rasqual/Y.tidy.bin \
	-k ../processed_data/rasqual/K.bin \
	-x ../processed_data/rasqual/X.bin \
	-n 52 -j $feat_id -l $n_rsnp -m $n_fsnp \
    -s $exon_start_positions -e $exon_end_positions \
    --imputation-quality 0.8 --imputation-quality-fsnp 0.8 \
    -f $gene_name -z --n_threads 1 \
    -c $tss -w $window_size -v \
    -r > ../processed_data/rasqual/output/${gene_id}_${gene_name}_w${window_size}_perm.txt
