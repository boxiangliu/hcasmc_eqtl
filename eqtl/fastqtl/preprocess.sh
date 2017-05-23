expr_dir=../processed_data/eqtl/fastqtl/expression/
cov_dir=../processed_data/eqtl/fastqtl/covariates/
if [[ ! -d $expr_dir ]]; then mkdir -p $expr_dir; fi
if [[ ! -d $cov_dir ]]; then mkdir -p $cov_dir; fi


# Prepare expression data for fastQTL:
Rscript eqtl/matrix_eqtl/gen_gene_loc.R ../../shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf $expr_dir/gene_loc.txt
Rscript eqtl/peer/normalize_rpkm.R ../processed_data/rnaseq/preprocess/combine_rpkm/combined.filter.rpkm $expr_dir/combined.filter.norm.rpkm
Rscript eqtl/fastqtl/gen_bed.R $expr_dir/gene_loc.txt $expr_dir/combined.filter.norm.rpkm $expr_dir/combined.filter.norm.bed
bgzip $expr_dir/combined.filter.norm.bed
tabix -p bed $expr_dir/combined.filter.norm.bed.gz


# Prepare covariates for fastQTL: 
Rscript eqtl/utils/combine_covariates.R \
--genotype_pc=../processed_data/genotype/genotype_pc/genotype_pcs.52samples.tsv \
--peer=../processed_data/eqtl/peer/factors.tsv \
--sample_info=../data/sample_info/sample_info.xlsx \
--output=$cov_dir/covariates.tsv \
--gender_coding=letter \
--num_geno_pc=4 \
--num_peer_factor=8

bgzip $cov_dir/covariates.tsv

