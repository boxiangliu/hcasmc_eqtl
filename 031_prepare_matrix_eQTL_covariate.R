#!/usr/bin/env Rscript
# bosh liu
# 2016/05/08
# durga
# prepare covariates data in matrix eQTL format

# library:
library(XLConnect)

# load covariates:
sample_sheet=readWorksheet(loadWorkbook("../processed_data/rna_wgs_match.reduced_050616.xlsx"),sheet=2)
sample_sheet=as.data.table(sample_sheet)


# cast RNA and DNA into strings: 
sample_sheet$RNA=as.character(sample_sheet$RNA)
sample_sheet$DNA=as.character(sample_sheet$DNA)


# remove 9052004_dase: 
sample_sheet=sample_sheet[RNA!='9052004_dase']


# sort sample_sheet alphanumerically: 
sample_sheet=sample_sheet[order(sample_sheet$DNA)]


# convert sample sheet to design matrix:
design=model.matrix(DNA~Sex+Manufacturer+Age_Imputed+AFR+EAS+AMR+EUR,sample_sheet)


# remove intercept: 
design=design[,colnames(design)!='(Intercept)']


# add rownames: 
rownames(design)=sample_sheet$DNA


# transpose: 
design=t(design)


# cast to data.table: 
design=as.data.table(design,keep.rownames=T)


# update the first column name: 
setnames(design,'rn','id')


# write table: 
write.table(design,'../processed_data/031_prepare_matrix_eQTL_covariate/covariates.txt',row.names=F,quote=F,sep='\t')