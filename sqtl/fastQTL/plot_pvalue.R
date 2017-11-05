#!/usr/bin/env Rscript
# bosh liu
# durga
# make qqplot of cis eqtl p-values

# library: 
library(data.table)
library(gap)
source('utils.R')
library(stringr)

# commandargs: 
args=commandArgs(T)
pval_file=args[1]
fig_dir=args[2]
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

# read input:
if (str_detect(pval_file,'.gz')) {
	pval=fread(sprintf('zcat %s',pval_file),header=F)
} else {
	pval=fread(pval_file,header=F)
}



# setnames: 
setnames(pval,c('intron','snp','distance','p-value','beta','varbeta'))


# randomly sample 100000 points:
set.seed(3)
pval_sample=pval[sample(nrow(pval),1000000),]


# make qqplot:
png(sprintf('%s/qqplot.png',fig_dir))
qqunif(pval_sample[,"p-value",with=F],pch=19)
dev.off()


# make histogram:
pdf(sprintf('%s/hist.pdf',fig_dir))
hist(unlist(pval_sample[,'p-value',with=F]),main='Nominal p-values',xlab='Nominal p-values',breaks=1000)
dev.off()
