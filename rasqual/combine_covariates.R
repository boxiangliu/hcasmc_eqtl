#!/usr/bin/env Rscript
# bosh liu
# durga
# combine gender, PEER factors and genotype PCs into one covariate file

# library:
library('R.utils')
library('XLConnect')
library(data.table)
library(dplyr)

# command arguments: 
args=commandArgs(trailingOnly=T,asValues=T,defaults=list(num_geno_pc=3,num_peer_factor=15))
genotype_pc_file=args$genotype_pc
peer_file=args$peer
sample_sheet_file=args$sample_info
output_file=args$output
gender_coding=args$gender_coding
num_geno_pc=args$num_geno_pc
num_peer_factor=args$num_peer_factor
row_and_colnames=as.logical(args$row_and_colnames)

# genotype_pc_file='../processed_data/160519_genotype_PCA/genotype_pcs.52samples.tsv'
# peer_file='../processed_data/160527/factors.tsv'
# sample_sheet_file='/srv/persistent/bliu2/HCASMC_eQTL/data/sample_info/sample_info.xlsx'
# output_file='../processed_data/160530/covariates.tsv'



# read input: 
message('INFO - reading input...')
genotype_pc=read.table(genotype_pc_file,header=T,row.names=1,check.names=F)
peer=read.table(peer_file,header=T,row.names=1,check.names=F)
sample_sheet=readWorksheet(loadWorkbook(sample_sheet_file),sheet=5)
gender=t(sample_sheet$Sex)


# recode gender (M=0, F=1):
if (gender_coding=='numerical'){
	gender=ifelse(gender=="M",0,1)
} else if (gender_coding=='letter'){
	# do nothing
} else {stop("gender coding should be either 'numerical' or 'letter'")}


# Converting gender to data.frame:
gender=data.frame(gender)
colnames(gender)=sample_sheet$DNA
if(nrow(gender)!=1) stop('gender should contain only 1 row!')


# check whether columns are in the same order:
message('INFO - checking column orders...')


if(!all(colnames(genotype_pc)==colnames(peer))){
	warning('Columns in genotype PC and PEER in different order.\nChanging the order of genotype PC to match PEER.')
	genotype_pc=genotype_pc[,match(colnames(peer),colnames(genotype_pc))]
}

if(!all(colnames(peer)==names(gender))){
	warning('Columns in gender and PEER in different order.\nChanging the order of gender to match PEER.')
	gender=gender[,match(colnames(peer),colnames(gender))]
}

stopifnot(all.equal(colnames(genotype_pc),colnames(peer)),all.equal(colnames(peer),names(gender)))


# combine covariates:
message('INFO - combining covariates...')
covariates=rbind(genotype_pc[1:num_geno_pc,],peer[1:num_peer_factor,],gender)


# set column names:
rownames(covariates)=c(paste0("C",seq(num_geno_pc)),paste0('InferredCov',seq(num_peer_factor)),'gender')
covariates=as.data.table(covariates,keep.rownames=T)
setnames(covariates,'rn','id')
# covariates[,ID:=NULL]


# transpose covariates:
# covariates=t(covariates)


# write output:
message('INFO - checking output...')
if (!dir.exists(dirname(output_file))){dir.create(dirname(output_file),recursive=TRUE)}
if (row_and_colnames){
	fwrite(covariates,file=output_file,sep='\t')
} else {
	fwrite(covariates,file=output_file,sep='\t')
}
