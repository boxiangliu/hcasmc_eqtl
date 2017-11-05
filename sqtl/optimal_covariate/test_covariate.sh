geno_pc=$1
splice_pc=$2
out_dir=$3
mkdir -p $out_dir 

echo INFO - genotype PC = $geno_pc
echo INFO - splicing PC = $splice_pc


echo INFO - making covariates
Rscript rasqual/combine_covariates.R \
	--genotype_pc=../processed_data/genotype/genotype_pc/genotype_pcs.52samples.tsv \
	--peer=../data/rnaseq2/leafcutter_wasp/cluster/sqtl_perind.counts.gz.PCs \
	--sample_info=/srv/persistent/bliu2/HCASMC_eQTL/data/sample_info/sample_info.xlsx \
	--output=$out_dir/covariates-${geno_pc}_geno_pc-${splice_pc}_expr_pc.tsv \
	--gender_coding=letter \
	--num_geno_pc=$geno_pc \
	--num_peer_factor=$splice_pc \
	--row_and_colnames=TRUE


echo INFO - running FastQTL
/users/zappala/software/fastqtl/bin/fastQTL \
--vcf ../data/joint3/asvcf/phased_and_imputed.chr22.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz \
--bed ../data/rnaseq2/leafcutter_wasp/cluster/sqtl_perind.counts.gz.qqnorm_chr22.gz \
--cov $out_dir/covariates-${geno_pc}_geno_pc-${splice_pc}_expr_pc.tsv \
--region chr22 \
--window 1e5 \
--out $out_dir/chr22.nominal.geno_pc-$geno_pc.splice_pc-$splice_pc.txt.gz 


echo INFO - counting significant loci
Rscript sqtl/optimal_covariate/calc_FDR.R \
$out_dir/chr22.nominal.geno_pc-$geno_pc.splice_pc-$splice_pc.txt.gz \
$out_dir/chr22.nominal.geno_pc-$geno_pc.splice_pc-$splice_pc.sig.txt \
$geno_pc $splice_pc
