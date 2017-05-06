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
setnames(fastqtl, c('gene_id','num_variants','shape1','shape2','dummy','best_variant','distance','nominal_p','slope','permute_p','beta_p'))


# perform bonferroni correction:
fastqtl[,bonferroni:=p.adjust(beta_p,method='bonferroni')]


# perform bh correction:
fastqtl[,bh:=p.adjust(beta_p,method='BH')]


# perform q-value correction: 
fastqtl$qvalues=qvalue(fastqtl$beta_p)$qvalues


# write output:
write.table(fastqtl,output_file,quote=F,sep='\t',row.names=F,col.names=T)

