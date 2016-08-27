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
# fastqtl_file='../processed_data/160805/hcasmc.eqtl.pc4.peer8.perm.txt'
# output_file='../processed_data/160805/hcasmc.eqtl.pc4.peer8.perm.padj.txt'

# read input: 
fastqtl=fread(fastqtl_file,header=F)


# column names:
# see code here: https://github.com/francois-a/fastqtl/blob/master/src/analysisPermutation.cpp line 97-139
# pheno
# nvar
# shape1
# shape2
# true df
# geno
# dist
# genotype_ma_samples
# genotype_ma_count
# maf
# genotype_ref_factor
# npval
# slope
# se
# ppval
# bpval
setnames(fastqtl, c("pid", "nvar", "shape1", "shape2", "df", "sid", "dist","ma_samples","ma_count","maf","ref","npval", "slope","se","ppval", "bpval"))

# perform bonferroni correction:
fastqtl[,bf:=p.adjust(bpval,method='bonferroni')]


# perform bh correction:
fastqtl[,bh:=p.adjust(bpval,method='BH')]


# perform q-value correction: 
fastqtl$qval=qvalue(fastqtl$bpval)$qvalues


# write output:
write.table(fastqtl,output_file,quote=F,sep='\t',row.names=F,col.names=T)

