#!/usr/bin/env Rscript
# bosh liu
# 2016/05/08
# durga
# generate snps location file

library(stringr)
library(data.table)

# command args: 
args=commandArgs(T)
input=args[1]
output=args[2]


# read genotype file:
genotypes=fread(input,header=T)


# extract genotype location:
split_id=str_split_fixed(genotypes$id,"_",4)


# construct genotype location file:
genotype_loc=data.frame(id=genotypes$id,chr=split_id[,1],pos=split_id[,2])


# write table: 
write.table(genotype_loc,file=output,row.names=F,quote=F,sep='\t')