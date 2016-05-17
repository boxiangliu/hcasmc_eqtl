#!/usr/bin/env Rscript
# bosh liu
# 2016/05/08
# durga
# keep genotypes with maf > 0.05

# command args: 
args=commandArgs(T)
input=args[1]
output=args[2]


# read genotype:
genotypes=fread(input,header=T)


# calculate allele frequencies: 
allele_count=rowSums(genotypes[,-c('id'),with=F])


# total number of alleles:
total_allele=(ncol(genotypes)-1)*2


# calculate minor allele frequency:
allele_freq=allele_count/total_allele
maf=pmin(allele_freq,1-allele_freq)


# subset genotype to maf > 0.05:
genotypes=genotypes[maf>=0.05]


# output table:
write.table(genotypes,file=output,row.names=F,quote=F,sep='\t')


