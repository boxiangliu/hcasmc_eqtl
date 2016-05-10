#!/usr/bin/env Rscript
# bosh liu
# 2016/05/08
# durga
# keep genotypes with maf > 0.05

library(data.table)
source('utils.R')

# command args: 
args=commandArgs(T)
chr=args[1]


# paths: 
output_dir='../processed_data/031_subset_genotype_by_maf'


# read genotype:
genotype_file=sprintf('../processed_data/031_prepare_matrix_eQTL_genotype/chr%s.genotype.txt',chr)
message(genotype_file)
genotypes=fread(genotype_file,header=T)


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
write.table(genotypes,file=sprintf('%s/chr%s.genotype.txt',output_dir,chr),row.names=F,quote=F,sep='\t')


