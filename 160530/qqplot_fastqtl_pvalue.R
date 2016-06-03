#!/usr/bin/env Rscript
# bosh liu
# durga
# plot nominal fastqtl p-values 


# library:
library(gap)
library(data.table)


# command args: 
args=commandArgs(T)
fastqtl_file=args[1]
hist_figure=args[2]
qqplot_figure=args[3]



# read input: 
fastqtl=fread(fastqtl_file,header=F)


# set column names:
setnames(fastqtl, c('gene_id','num_variants','shape1','shape2','dummy','best_variant','distance','nominal_p','slope','permute_p','beta_p'))


# make histogram of p-value:
pdf(hist_figure)
hist(fastqtl$nominal_p,main='Nominal p-value',xlab='p-value')
hist(fastqtl$permute_p,main='P-value via direct permutation',xlab='p-value')
hist(fastqtl$beta_p,main='P-value via beta approximation',xlab='p-value')
dev.off()


# make qqplot of p-value: 
pdf(qqplot_figure)
qqunif(fastqtl$nominal_p,main='nominal p-value')
qqunif(fastqtl$permute_p,main='P-value via direct permutation')
qqunif(fastqtl$beta_p,main='P-value via beta approximation')
dev.off()




