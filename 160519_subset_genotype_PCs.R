#!/bin/R
# Bosh Liu
# 2016/05/19
# durga
# subset genotype PCs to individuals with RNAseq

# library:
library('XLConnect')


# command args: 
args=commandArgs(T)
input_file=args[1]
input_file='../processed_data/160519_genotype_PCA/genotype_pcs.tsv'
sample_sheet_file=args[2]
sample_sheet_file='/srv/persistent/bliu2/HCASMC_eQTL/data/sample_info/sample_info.xlsx'
output_file=args[3]
output_file='../processed_data/160519_genotype_PCA/genotype_pcs.52samples.tsv'

# read input:
input=read.table(input_file,header=T,row.names=1,check.names=F)
sample_sheet=readWorksheet(loadWorkbook(sample_sheet_file),sheet=5)


# read input:
samples=sample_sheet$DNA
idx=match(samples,colnames(input))
output=input[,idx]


# write output: 
write.table(output,output_file,quote=F,sep='\t',row.names=T,col.names=T)