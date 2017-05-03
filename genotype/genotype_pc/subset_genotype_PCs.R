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
sample_sheet_file=args[2]
output_file=args[3]

# read input:
input=read.table(input_file,header=T,row.names=1,check.names=F)
sample_sheet=readWorksheet(loadWorkbook(sample_sheet_file),sheet=5)


# read input:
samples=sample_sheet$DNA
idx=match(samples,colnames(input))
output=input[,idx]


# write output: 
write.table(output,output_file,quote=F,sep='\t',row.names=T,col.names=T)