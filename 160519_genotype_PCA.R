#!/bin/R
#r v 3.2.3
# Brunilda Balliu
# modified by Bosh Liu
# 2016/05/19
# durga
# get genotype PCs

# paths:
figure_path='../figures/160519_genotype_PCA/'

# source("http://bioconductor.org/biocLite.R")
# biocLite("gdsfmt")
# biocLite("SNPRelate")
library('SNPRelate')
library('XLConnect')

# Convert the PLINK files to the GDS file
# Data should be QC-ed for missing rate, MAF and HWE!
# bed.fn <- "../processed_data/160519_genotype_PCA/recalibrated_biallelic_SNP.beagle.rename.bed"
# fam.fn <- "../processed_data/160519_genotype_PCA/recalibrated_biallelic_SNP.beagle.rename.fam"
# bim.fn <- "../processed_data/160519_genotype_PCA/recalibrated_biallelic_SNP.beagle.rename.bim"
bed.fn <- "../processed_data/160519_genotype_PCA/recalibrated_biallelic_SNP.beagle.rename.dr2.bed"
fam.fn <- "../processed_data/160519_genotype_PCA/recalibrated_biallelic_SNP.beagle.rename.dr2.fam"
bim.fn <- "../processed_data/160519_genotype_PCA/recalibrated_biallelic_SNP.beagle.rename.dr2.bim"

# convert
# snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "../processed_data/160519_genotype_PCA/recalibrated_biallelic_SNP.beagle.rename.gds")
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "../processed_data/160519_genotype_PCA/recalibrated_biallelic_SNP.beagle.rename.dr2.gds")
# genofile <- snpgdsOpen("../processed_data/160519_genotype_PCA/recalibrated_biallelic_SNP.beagle.rename.gds")
genofile <- snpgdsOpen("../processed_data/160519_genotype_PCA/recalibrated_biallelic_SNP.beagle.rename.dr2.gds")
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, maf=0.05)
snpset.id <- unlist(snpset)

# Run pca
pca <- snpgdsPCA(genofile, snp.id=snpset.id, maf=0.05, num.thread=2)



# read genomic ethnicity: 
sample_sheet=readWorksheet(loadWorkbook("../processed_data/rna_wgs_match.reduced_050616.xlsx"),sheet=4)
sample_sheet=sample_sheet[match(pca$sample.id,sample_sheet$DNA),]


#Plot PCA:
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE,
                  pop = factor(sample_sheet$Genomic_Ethnicity))
pdf(paste0(figure_path,'genotype_PC1_and_PC2.pdf'))
plot(tab$EV1, tab$EV2, col= as.integer(tab$pop), xlab="eigenvector 1", ylab="eigenvector 2")
legend("bottomright", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))
dev.off()


# variance proportion:
pdf(paste0(figure_path,'genotype_PC_pairs.pdf'),width=12,height=12)
pc.percent <- pca$varprop*100
lbls <- paste("PC", 1:5, "\n", format(pc.percent[1:5], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:5], labels=lbls, col=as.integer(tab$pop))
dev.off()

