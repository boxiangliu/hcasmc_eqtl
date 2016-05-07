#!/usr/bin/env Rscript
# bosh liu
# 2016/05/07
# durga
# create symbolic link from old RNAseq name to new RNAseq names
# using ../processed_data/rna_wgs_match.reduced_050616.xlsx

# load libraries
library(XLConnect)


# read sample sheet
sample_sheet=readWorksheet(loadWorkbook("../processed_data/rna_wgs_match.reduced_050616.xlsx"),sheet=1)


# keep only columns for old and new RNAseq sample names:
sample_sheet=sample_sheet[,c('rna','RNA.New.Name')]


# rename columns for easy of reference: 
colnames(sample_sheet) = c('rna_old','rna_new')


# convert to sample_sheet data.table for easy data manipulation:
sample_sheet=as.data.table(sample_sheet)


# subset to rows whose rna_new column entries are not NA:
sample_sheet=sample_sheet[!is.na(rna_new)]


# create symbolic link from old to new RNAseq name 
# using reverse strand counts because 
# only 90702_Nextseq was sequenced on the forward strand
# but it will not be used for eQTL mapping
src_dir=normalizePath('../data/rnaseq/expression/reverse/')
dst_dir=normalizePath('../data/rnaseq/expression/working_set/')
if (!dir.exists(dst_dir)) {dir.create(dst_dir)}
for (i in seq(nrow(sample_sheet))) {
	rna_old=sample_sheet[i,rna_old]
	rna_new=sample_sheet[i,rna_new]
	message('create symbolic link from ',rna_old,' to ',rna_new)
	file.symlink(from=paste(src_dir,rna_old,sep='/'),to=paste(dst_dir,rna_new,sep='/'))
}

