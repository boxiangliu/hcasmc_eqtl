#!/usr/bin/env Rscript 
# bosh liu
# apply bonferroni and BH adjustment procedure:

# library:
library(data.table)
library(qvalue)


# command args: 
args=commandArgs(T)
fastqtl_file=args[1]
output_file=args[2]


# read input: 
fastqtl=fread(fastqtl_file,header=F)


# setnames:
setnames(fastqtl, c('gene_id','sid','distance','pval','beta','varbeta'))


# perform bh correction:
fastqtl[,padj:=p.adjust(pval,method='BH')]


# perform q-value correction: 
fastqtl$qval=qvalue(fastqtl$pval)$qvalues


# write output:
write.table(fastqtl,output_file,quote=F,sep='\t',row.names=F,col.names=T)

