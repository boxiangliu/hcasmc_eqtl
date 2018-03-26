# Genotype PCs: 
# convert VCF to plink BED file:
plink --vcf ../data/joint2/recalibrated_biallelic_SNP.beagle.rename.dr2.vcf --keep-allele-order --make-bed --out ../processed_data/genotype/genotype_pc/recalibrated_biallelic_SNP.beagle.rename.dr2

# find genotype PCs: 
Rscript genotype/genotype_pc/genotype_PCA.R

# create sample_info directory and sample_list.txt:
mkdir /srv/persistent/bliu2/HCASMC_eQTL/data/sample_info
ln /srv/persistent/bliu2/HCASMC_eQTL/processed_data/rna_wgs_match.reduced_050616.xlsx /srv/persistent/bliu2/HCASMC_eQTL/data/sample_info/sample_info.xlsx

# subset to 52 individuals with RNAseq sample:
Rscript genotype/genotype_pc/subset_genotype_PCs.R ../processed_data/genotype/genotype_pc/genotype_pcs.tsv ../data/sample_info/sample_info.xlsx ../processed_data/genotype/genotype_pc/genotype_pcs.52samples.tsv
Rscript genotype/genotype_pc/make_PCA_plot.R 

# Comment:
# seems that 1020301 is an outlier? It has low heterozygosity rate (from the plink/seq analysis)
