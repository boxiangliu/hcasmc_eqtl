# setup:
scripts=eqtl/matrix_eqtl
processed_data=../processed_data/160530/
figures=../figures/160530/
mkdir $scripts $processed_data $figures

# prepare genotype data
# filter minor allele frequency 0.05: 
bcftools view \
	--min-af 0.05 --max-af 0.95 \
	-Ov -o ../data/joint2/recalibrated_biallelic_SNP.beagle.rename.dr2.hwe.maf.vcf \
	../data/joint2/recalibrated_biallelic_SNP.beagle.rename.dr2.hwe.vcf

bcftools query -H \
	-S ../data/joint2/sample_list.txt \
	-f '%CHROM\_%POS\_%REF\_%ALT[\t%DS]\n' \
	-o $processed_data/dosage.tsv \
	../data/joint2/recalibrated_biallelic_SNP.beagle.rename.dr2.hwe.maf.vcf

cp $processed_data/dosage.tsv $processed_data/dosage.tmp

sed -e "s/# \[1\]CHROM_\[2\]POS_\[3\]REF_\[4\]ALT/id/" \
	-e "s/\[[[:digit:]]\+\]//g" -e "s/:DS//g" \
	$processed_data/dosage.tmp > $processed_data/dosage.tsv
rm $processed_data/dosage.tmp


# prepare expression data: 
ln ../processed_data/160527/combined.filter.rpkm $processed_data/combined.filter.rpkm
Rscript 160527/normalize_rpkm.R $processed_data/combined.filter.rpkm $processed_data/combined.filter.norm.rpkm


# prepare genotype location: 
cp 031_gen_snps_loc.R $scripts/gen_snps_loc.R
Rscript $scripts/gen_snps_loc.R $processed_data/dosage.tsv $processed_data/genotype_loc.txt


# prepare gene location: 
cp 031_gen_gene_loc.R $scripts/gen_gene_loc.R
subl $scripts/gen_gene_loc.R
Rscript $scripts/gen_gene_loc.R /srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf $processed_data/gene_loc.txt


# find optimum number of PEER factors:
mkdir $processed_data/find_optimum_num_PEER_factors_matrixeqtl
bash $scripts/find_optimal_num_PEER_factors_matrix_eQTL.sh 
bash $scripts/run_count_num_sig_association.sh 
Rscript $scripts/plot_num_eqtl_vs_cov.R \
	$processed_data/find_optimum_num_PEER_factors_matrixeqtl/num_eqtls_vs_cov.fdr.txt \
	$figures/num_eqtls_vs_cov.fdr.pdf


# run matrix eQTL:
cp 160516_matrix_eQTL.R $scripts/run_matrix_eQTL.R
# Rscript $scripts/run_matrix_eQTL.R \
# 	$processed_data/dosage.tsv \
# 	$processed_data/genotype_loc.txt \
# 	$processed_data/combined.filter.norm.rpkm \
# 	$processed_data/gene_loc.txt \
# 	$processed_data/covariates.tsv \
# 	$processed_data/
# Use $processed_data/find_optimum_num_PEER_factors_matrixeqtl/pc3.peer8.2.cis.txt


# make qqplot and histogram of p-value distribution: 
mkdir ../figures/160530
Rscript $scripts/qqplot_pvalue.R ../processed_data/160530/cis.txt ../figures/160530/

