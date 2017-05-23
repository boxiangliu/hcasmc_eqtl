vcf_dir=../data/joint3/asvcf_sid/
expr_dir=../processed_data/eqtl/fastqtl/expression/
cov_dir=../processed_data/eqtl/fastqtl/covariates/
out_dir=../processed_data/eqtl/fastqtl/output/nominal/
log_dir=../logs/fasqtl/nominal/
if [[ ! -d $out_dir ]]; then mkdir -p $out_dir; fi
if [[ ! -d $log_dir ]]; then mkdir -p $log_dir; fi


parallel -j10 ../../tools/fastqtl/bin/fastQTL \
	--vcf $vcf_dir/phased_and_imputed.chr{}.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz \
	--bed $expr_dir/combined.filter.norm.bed.gz \
	--out $out_dir/chr{}.txt.gz \
	--cov $cov_dir/covariates.tsv.gz \
	--normal \
	--maf-threshold 0.05 \
	--region chr{} '>' $log_dir/chr{}.log ::: `seq 22`

zcat $out_dir/chr{1..22}.txt.gz | bgzip > $out_dir/all.txt.gz
