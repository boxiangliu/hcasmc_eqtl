#!/usr/bin/env Rscript
# bosh liu
# 2016/05/06
# durga
# compare the imputed genotype with original genotypes

impute2_filename='../data/joint/imputation/chr22.genotype'
imputed=fread(impute2_filename)
# imputed=imputed[,V1:=NULL]
# setnames(imputed,'V2','SNP')

# read original genotype:
original_filename='../processed_data/028_imputed_genotype/recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.chr22.vcf'
original=fread(original_filename)
setnames(original,'#CHROM','CHROM')
# setnames(imputed,c('SNP',colnames(original)[10:76]))

# compare dimension:
nrow(impute2) # 1110217
nrow(original) # 158561


# what percentage of original genotypes are not in imputed genotypes? [90%]
imputed_SNP=imputed[,SNP]
original_SNP=original[,paste(CHROM,POS,REF,ALT,sep="_")]
original$SNP=original_SNP
temp=original_SNP %in% imputed_SNP
mean(temp) # 0.9009214


# are imputed variants same as the original variants? [YES]
head(temp) # [1] FALSE  TRUE FALSE FALSE  TRUE  TRUE
imputed[SNP==original_SNP[2]]
original[SNP==original_SNP[2]]

