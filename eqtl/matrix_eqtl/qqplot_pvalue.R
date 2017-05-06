#!/usr/bin/env Rscript
# bosh liu
# durga
# make qqplot of cis eqtl p-values

# library: 
library(gap)
source('utils.R')


# commandargs: 
args=commandArgs(T)
pval_file=args[1]
figure_prefix=args[2]
# pval_file='../processed_data/160530/cis.txt'
# figure_prefix='../figures/160530/'


# read input:
pval=fread(pval_file,header=T)


# randomly sample 100000 points:
set.seed(3)
pval_sample=pval[sample(nrow(pval),100000),]


# make qqplot:
pdf(paste0(figure_prefix,'qqplot.pdf'))
qqunif(pval_sample[,"p-value",with=F],pch=19)
dev.off()


# make histogram:
pdf(paste0(figure_prefix,'histogram.pdf'))
hist(unlist(pval[,'p-value',with=F]),main='Nominal p-values',xlab='Nominal p-values')
dev.off()
