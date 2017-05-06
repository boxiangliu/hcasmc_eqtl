scripts=eqtl/fastqtl/
processed_data=../processed_data/160530/
figures=../figures/160530

# Prepare genotype data for fastQTL: 
bcftools annotate -Oz -o /srv/persistent/bliu2/HCASMC_eQTL/data/joint2/recalibrated_biallelic_SNP.beagle.rename.dr2.hwe.maf.vcf.id.gz --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' /srv/persistent/bliu2/HCASMC_eQTL/data/joint2/recalibrated_biallelic_SNP.beagle.rename.dr2.hwe.maf.vcf.gz
tabix -p vcf /srv/persistent/bliu2/HCASMC_eQTL/data/joint2/recalibrated_biallelic_SNP.beagle.rename.dr2.hwe.maf.vcf.id.gz


# Prepare expression data for fastQTL: 
Rscript $scripts/gen_bed.R $processed_data/gene_loc.txt $processed_data/combined.filter.norm.rpkm $processed_data/combined.filter.norm.bed
bgzip $processed_data/combined.filter.norm.bed
tabix -p bed $processed_data/combined.filter.norm.bed.gz


# Prepare covariates for fastQTL: 
Rscript eqtl/utils/combine_covariates.R --genotype_pc=../processed_data/160519_genotype_PCA/genotype_pcs.52samples.tsv --peer=../processed_data/160527/factors.tsv --sample_info=/srv/persistent/bliu2/HCASMC_eQTL/data/sample_info/sample_info.xlsx --output=$processed_data/covariates.fastqtl.tsv --gender_coding=letter
bgzip $processed_data/covariates.fastqtl.tsv


# Run fastQTL (nominal pass):
n=0
for i in {1..22}; do 
n=$(($n+1))
if [[ n -gt 10 ]]; then wait; n=0; fi 
	bash $scripts/run_fastqtl.nominal.sh /srv/persistent/bliu2/HCASMC_eQTL/data/joint2/recalibrated_biallelic_SNP.beagle.rename.dr2.hwe.maf.vcf.id.gz $processed_data/combined.filter.norm.bed.gz $processed_data/find_optimum_num_PEER_factors/covariates.fastqtl.pc3.peer8.tsv.gz $processed_data/fastqtl_nominal/fastqtl.chr$i.pc3.peer8.txt.gz chr$i > $processed_data/fastqtl_nominal/run_fastqtl.chr$i.log &
done 
zcat $processed_data/fastqtl_nominal/fastqtl.chr{1..22}.pc3.peer8.txt.gz > $processed_data/fastqtl_nominal/fastqtl.allpairs.pc3.peer8.txt


# Create symbolic link from fastQTL result of the nominal pass to the data directory: 
ln -s /srv/persistent/bliu2/HCASMC_eQTL/processed_data/fastqtl_nominal/ /srv/persistent/bliu2/HCASMC_eQTL/data/eQTL/fastqtl/

# Adjust p-value for fastQTL result from nominal pass: 
Rscript $scripts/fastqtl_pvalue_corrections.nominal_pass.R $data/eQTL/fastqtl/fastqtl_nominal/fastqtl.allpairs.pc3.peer8.txt $data/eQTL/fastqtl/fastqtl_nominal/fastqtl.allpairs.pc3.peer8.padj.txt


# Make histogram and qqplot of fastQTL nominal p-values:
Rscript $scripts/qqplot_fastqtl_pvalue.R $processed_data/fastqtl.txt $figures/fastqtl_histogram.pdf $figures/fastqtl_qqplot.pdf


# find optimal number of PEER factors:
mkdir $processed_data/find_optimum_num_PEER_factors/
bash $scripts/find_optimal_num_PEER_factors.sh
Rscript $scripts/plot_num_egene_vs_cov.R $figures/num_egene_vs_cov.pdf


# Copy the optimal eQTL set and remove the rest: 
mkdir $processed_data/fastqtl_perm/
cp $processed_data/find_optimum_num_PEER_factors/fastqtl.pc4.peer8* $processed_data/fastqtl_perm/
rm -r $processed_data/find_optimum_num_PEER_factors/


# Create symbolic link from fastQTL result of the permutation pass to the data directory:
ln -s /srv/persistent/bliu2/HCASMC_eQTL/processed_data/fastqtl_perm/ /srv/persistent/bliu2/HCASMC_eQTL/data/eQTL/fastqtl/fastqtl_perm


# How does our discover compare with GTEx: 
Rscript $scripts/plot_num_egenes_vs_sample_size.R \
	$processed_data/gtex.v6p.egenes.summary.txt \
	$processed_data/find_optimum_num_PEER_factors/fastqtl.pc4.peer8.padj.txt \
	$figures/num_egenes_vs_sample_size.pdf


# num eGenes discovered by FDR: 
Rscript $scripts/plot_num_egene_by_fdr.R $processed_data/fastqtl.padj.txt $figures/num_egenes_by_fdr.pdf 


# investigate the top 5 eGenes:
Rscript $scripts/plot_expression_vs_genotype.R \
	$processed_data/combined.filter.norm.rpkm \
	$processed_data/dosage.tsv \
	$processed_data/covariates.tsv \
	$processed_data/fastqtl.padj.txt \
	$figures