#!/usr/bin/env Rscript 
# bosh liu
# apply bonferroni and BH adjustment procedure:

# library:
library(data.table)
library(qvalue)


# command args: 
args=commandArgs(T)
in_file=args[1]
output_file=args[2]


# read input: 
eqtl=fread(in_file,header=T)


# perform bh correction:
eqtl[,padj:=p.adjust(pval,method='BH')]


# perform q-value correction: 
eqtl$qval=qvalue(eqtl$pval)$qvalues


# write output:
write.table(eqtl,output_file,quote=F,sep='\t',row.names=F,col.names=T)

