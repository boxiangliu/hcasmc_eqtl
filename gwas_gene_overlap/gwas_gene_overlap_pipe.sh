Rscript gwas_gene_overlap/ldscore_regression/tissue_specific_gene.R
Rscript gwas_gene_overlap/ldscore_regression/tissue_specific_snp_annotation.R
bash gwas_gene_overlap/ldscore_regression/gwas_sumstats.sh

# bash gwas_gene_overlap/ldscore_regression/ldscore.sh
bash gwas_gene_overlap/ldscore_regression/ldscore_merged.sh

# bash gwas_gene_overlap/ldscore_regression/partition_heritability.sh
bash gwas_gene_overlap/ldscore_regression/partition_heritability_merged.sh \
../processed_data/gwas_gene_overlap/ldscore_regression/gwas_sumstats/cad.sumstats.gz \
../processed_data/gwas_gene_overlap/ldscore_regression/partition_heritability_merged/

bash gwas_gene_overlap/ldscore_regression/partition_heritability_merged.sh \
../processed_data/gwas_gene_overlap/ldscore_regression/gwas_sumstats/scz.sumstats.gz \
../processed_data/gwas_gene_overlap/ldscore_regression/partition_heritability_merged_scz/

Rscript gwas_gene_overlap/ldscore_regression/plot_scz.R