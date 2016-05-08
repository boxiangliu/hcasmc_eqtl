#!/usr/bin/env Rscript 
# bosh liu
# 2016/04/26
# durga
# generate a table of covariates 

# library:
library(XLConnect)


# samples: 
sample_sheet=readWorksheet(loadWorkbook("../processed_data/rna_wgs_match.reduced_050616.xlsx"),sheet=1)


# sex:
sex=read.table('../processed_data/sex.tsv',header=T)
temp=merge(sample_sheet,sex,by.x='RNA.New.Name',by.y='sample',all=T)
temp[!is.na(temp$RNA.New.Name),c('dna','RNA.New.Name', 'chrY','chrX','reported','concordant')]
# I manually entered sex into rna_wgs_match.reduced_050616.xlsx

# age: 
sample_sheet2=readWorksheet(loadWorkbook("../raw_data/R21R33\ Cell\ line\ summary.xlsx"),sheet=1)
sample_sheet2
temp=merge(sample_sheet,sample_sheet2,by.x='RNA.New.Name',by.y='Cell.Line')
temp=as.data.table(temp)
temp[,.(RNA.New.Name,Manufacturer,Catalog..,Age,Gender,Ethnicity)]
# I manually entered age into rna_wgs_match.reduced_050616.xlsx


# previously a separate sheet was made with only samples in the working set:
# read the new sheet: 
sample_sheet=readWorksheet(loadWorkbook("../processed_data/rna_wgs_match.reduced_050616.xlsx"),sheet=2)



# ancestry proportions:
ancestry_filename='../processed_data/ancestry_proportions.tsv'
ancestry=fread(ancestry_filename,header=T)
temp=merge(sample_sheet,ancestry,by.x='DNA',by.y='ID',all.x=T)
temp=as.data.table(temp)


# rows are ordered when merging
# revert row order to original order  
original_order=sample_sheet$DNA
idx=match(original_order,temp$DNA)
temp=temp[idx,]
temp[,.(AFR,EAS,AMR,EUR)]
# copied and pasted from R to excel
# use text import wizard (option: fixed width text)
# manually filled in the NA columns (1401 and 1508)
