#!/usr/bin/env Rscript
# bosh liu
# 2016/05/08
# durga
# prepare expression data in matrix eQTL format
# load expression data:

library(DESeq2)
# path: 
output_dir='../processed_data/031_prepare_matrix_eQTL_expression'


# read variance stabilized counts:
vsd=readRDS('../processed_data/030_variance_stabilize/vsd.rds')


# get expression matrix: 
expression=assay(vsd)


# remove "_dase" affix: 
colnames(expression)=str_split_fixed(colnames(expression),"_",n=2)[,1]



# check expression and snp are in sync:
stopifnot(colnames(expression)==sort(colnames(expression)))


# cast expression from a matrix to a data.table: 
expression=as.data.table(expression,keep.rownames=T)
setnames(expression,'rn','id')


# save tables:
if (!dir.exists(output_dir)){dir.create(output_dir)}
write.table(expression,file=paste(output_dir,'expression.txt',sep='/'),row.names=F,quote=F,sep='\t')
