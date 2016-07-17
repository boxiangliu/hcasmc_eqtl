scripts=./160530
processed_data=../processed_data/160530
num_geno_pc=3
num_peer_factor=15

covariates=$processed_data/find_optimum_num_PEER_factors_matrixeqtl/covariates.matrixeqtl.pc$num_geno_pc.peer$num_peer_factor.tsv
grep -e "chr22" -e "id" $processed_data/dosage.tsv > $processed_data/dosage.chr22.tsv
grep -e "chr22" -e "id" $processed_data/gene_loc.txt > $processed_data/gene_loc.chr22.txt
grep -e "chr22" -e "id" $processed_data/genotype_loc.txt > $processed_data/genotype_loc.chr22.txt
Rscript $scripts/run_matrix_eQTL.R \
$processed_data/dosage.chr22.tsv \
$processed_data/genotype_loc.chr22.txt \
$processed_data/combined.filter.norm.rpkm \
$processed_data/gene_loc.chr22.txt \
$covariates \
$processed_data/find_optimum_num_PEER_factors_matrixeqtl/pc$num_geno_pc.peer$num_peer_factor.2.