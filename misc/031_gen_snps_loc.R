#!/usr/bin/env Rscript
# bosh liu
# 2016/05/08
# durga
# generate snps location file


library(data.table)
source('utils.R')

# command args: 
args=commandArgs(T)
chr=args[1]


# paths: 
output_dir='../processed_data/031_gen_snps_loc/'


# read genotype file:
genotype_file=sprintf('../processed_data/031_subset_genotype_by_maf/chr%s.genotype.txt',chr)
message(genotype_file)
genotypes=fread(genotype_file,header=T)


# extract genotype location:
split_id=str_split_fixed(genotypes$id,"_",4)


# construct genotype location file:
genotype_loc=data.frame(id=genotypes$id,chr=split_id[,1],pos=split_id[,2])


# write table: 
write.table(genotype_loc,file=sprintf('%s/chr%s.genotype_loc.txt',output_dir,chr),row.names=F,quote=F,sep='\t')