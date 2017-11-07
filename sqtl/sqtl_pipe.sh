# Create directories:
mkdir -p ../processed_data/sqtl/ ../figures/sqtl/


# Convert BAM to junction files:
bash sqtl/leafcutter/run_bam2junc.sh


# Intron clustering:
python /srv/persistent/bliu2/tools/leafcutter/clustering/leafcutter_cluster.py \
	-j ../data/rnaseq2/wasp/juncfiles.txt \
	-r ../data/rnaseq2/leafcutter_wasp/cluster/ \
	-o sqtl


# Prepare FastQTL input:
python /srv/persistent/bliu2/tools/leafcutter/scripts/prepare_phenotype_table.py \
	../data/rnaseq2/leafcutter_wasp/cluster/sqtl_perind.counts.gz \
	-p 15


# Change GRCh37 to hg19 coordinates:
parallel -j1 \
cat ../data/rnaseq2/leafcutter_wasp/cluster/sqtl_perind.counts.gz.qqnorm_chr{} "|" \
python sqtl/utils/b37_to_hg19.py ">" ../data/rnaseq2/leafcutter_wasp/cluster/sqtl_perind.counts.gz.qqnorm_chr{}_2 ::: {1..22}

parallel -j1 \
mv ../data/rnaseq2/leafcutter_wasp/cluster/sqtl_perind.counts.gz.qqnorm_chr{}_2 \
../data/rnaseq2/leafcutter_wasp/cluster/sqtl_perind.counts.gz.qqnorm_chr{} ::: {1..22}


# Bgzip and index the bed files:
bash ../data/rnaseq2/leafcutter_wasp/cluster/sqtl_perind.counts.gz_prepare.sh # bgzip and index


# Select optimal combination of covariates:
parallel -j10 --joblog ../processed_data/sqtl/optimal_covariate/test_covariate/log \
bash sqtl/optimal_covariate/test_covariate.sh \
{1} {2} ../processed_data/sqtl/optimal_covariate/test_covariate/ \
::: 2 3 4 ::: 1 2 3 4 5 6 7 8 9 12 15

cat ../processed_data/sqtl/optimal_covariate/test_covariate/chr22.nominal.*.sig.txt \
> ../processed_data/sqtl/optimal_covariate/test_covariate/sig.txt

Rscript sqtl/optimal_covariate/compare_covariate.R

# Add ancestry PC and sex:
Rscript rasqual/combine_covariates.R \
	--genotype_pc=../processed_data/genotype/genotype_pc/genotype_pcs.52samples.tsv \
	--peer=../data/rnaseq2/leafcutter_wasp/cluster/sqtl_perind.counts.gz.PCs \
	--sample_info=/srv/persistent/bliu2/HCASMC_eQTL/data/sample_info/sample_info.xlsx \
	--output=../processed_data/sqtl/covariate/covariates-3_geno_pc-6_splice_pc.tsv \
	--gender_coding=letter \
	--num_geno_pc=3 \
	--num_peer_factor=6 \
	--row_and_colnames=TRUE


# sQTL mapping with fastQTL in nominal mode:
parallel -j10 /users/zappala/software/fastqtl/bin/fastQTL \
--vcf ../data/joint3/asvcf/phased_and_imputed.chr{}.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz \
--bed ../data/rnaseq2/leafcutter_wasp/cluster/sqtl_perind.counts.gz.qqnorm_chr{}.gz \
--cov ../processed_data/sqtl/covariate/covariates-3_geno_pc-6_splice_pc.tsv \
--region chr{} \
--window 1e5 \
--out ../processed_data/sqtl/fastQTL/nominal/chr{}.nominal.txt.gz ::: {1..22}

zcat ../processed_data/sqtl/fastQTL/nominal/chr{1..22}.nominal.txt.gz | \
gzip > ../processed_data/sqtl/fastQTL/nominal/all.nominal.txt.gz

# sQTL mapping with fastQTL in permutation mode:
parallel -j10 /users/zappala/software/fastqtl/bin/fastQTL \
--vcf ../data/joint3/asvcf/phased_and_imputed.chr{}.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz \
--bed ../data/rnaseq2/leafcutter_wasp/cluster/sqtl_perind.counts.gz.qqnorm_chr{}.gz \
--cov ../processed_data/sqtl/covariate/covariates-3_geno_pc-6_splice_pc.tsv \
--region chr{} \
--window 1e5 \
--permute 1000 100000 \
--out ../processed_data/sqtl/fastQTL/permutation/chr{}.permutation.txt.gz ::: 1
{1..22}

zcat ../processed_data/sqtl/fastQTL/permutation/chr{1..22}.permutation.txt.gz | \
gzip > ../processed_data/sqtl/fastQTL/permutation/all.permutation.txt.gz

# Plot p-values: 
Rscript sqtl/fastQTL/plot_pvalue.R \
../processed_data/sqtl/fastQTL/nominal/all.nominal.txt.gz \
../figures/sqtl/fastQTL/plot_pvalue/


# Plot sQTL vs distance to TSS:
Rscript sqtl/fastQTL/plot_sqtl_vs_distance.R \
../processed_data/sqtl/fastQTL/nominal/all.nominal.txt.gz \
../figures/sqtl/fastQTL/plot_sqtl_vs_distance/


# Count sQTLs: 
Rscript sqtl/fastQTL/count_sqtl.R \
../processed_data/sqtl/fastQTL/nominal/all.nominal.txt.gz \
../processed_data/sqtl/fastQTL/count_sqtl/ \
fastqtl.nominal

Rscript sqtl/fastQTL/count_sqtl.R \
../processed_data/sqtl/fastQTL/permutation/all.permutation.txt.gz \
../processed_data/sqtl/fastQTL/count_sqtl/ \
fastqtl.permutation
