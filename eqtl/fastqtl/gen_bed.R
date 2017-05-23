#!/usr/bin/env Rscript
# bosh liu
# generate input bed file for fastQTL
library(data.table)

# command args: 
args=commandArgs(T)
gene_loc_file=args[1]
expression_file=args[2]
output_file=args[3]

# read inputs:  
gene_loc=fread(gene_loc_file,header=T)
expression=fread(expression_file,header=T)


# merge gene location and expression:
idx=match(expression$Name,gene_loc$id)
gene_loc_short=gene_loc[idx,]
stopifnot(gene_loc_short$id==expression$Name)
bed=data.table(gene_loc_short[,.(chr,left,right)],expression)


# set column names:
setnames(bed,c('chr', 'left','right','Name'),c('#chr', 'start','end','ID'))


# write output:
write.table(bed,output_file,quote=F,sep='\t',row.names=F,col.names=T)