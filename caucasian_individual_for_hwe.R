#!/usr/bin/env Rscript
# bosh liu
# 2016/05/15
# durga
# make sample file for Caucasian individuals for hwe filtering: 


# library:
library(XLConnect)


# paths:
input_file='../processed_data/rna_wgs_match.reduced_050616.xlsx'
output_file='../data/joint2/caucasian_for_hwe.txt'


# read input:
input=readWorksheet(loadWorkbook(input_file),sheet=4)
input

# select caucasian samples:
sample_names=input[input$Genomic_Ethnicity=='Caucasian','DNA']
sample_names=sort(as.character(sample_names))


# write sample names to table:
write.table(sample_names,file=output_file,quote=F,row.names=F,col.names=F)
