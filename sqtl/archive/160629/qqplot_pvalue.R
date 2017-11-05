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
qqplot_figure=args[2]
histogram_figure=args[3]
# pval_file='../processed_data/160629/sqtl.nominal.allpairs.normal.1e5.txt'
# qqplot_figure='../figures/160629/qqplot.normal.1e5.pdf'
# histogram_figure='../figures/160629/histogram.normal.1e5.pdf'


# read input:
pval=fread(pval_file,header=F)


# setnames: 
setnames(pval,c('intron','snp','distance','p-value','beta','varbeta'))


# randomly sample 100000 points:
set.seed(3)
pval_sample=pval[sample(nrow(pval),1000000),]


# make qqplot:
png(qqplot_figure)
qqunif(pval_sample[,"p-value",with=F],pch=19)
dev.off()


# make histogram:
pdf(histogram_figure)
hist(unlist(pval_sample[,'p-value',with=F]),main='Nominal p-values',xlab='Nominal p-values',breaks=1000)
dev.off()
