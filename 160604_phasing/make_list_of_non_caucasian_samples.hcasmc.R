#!/usr/bin/env Rscript
# bosh liu
# 2016/05/15
# durga
# make sample file for individuals in each ethnicity 


# library:
library(XLConnect)


# command args:
args=commandArgs(T)
input_file=args[1]
output_file=args[2]
# input_file='../processed_data/rna_wgs_match.reduced_050616.xlsx'
# output_file='../processed_data/160526/detect_sample_contamination_model_based'


# read input:
input=readWorksheet(loadWorkbook(input_file),sheet=4)


# select ethnic samples:
ethnicity=c('AA','Hispanic','Asian')
sample_names=input[input$Genomic_Ethnicity!='Caucasian','DNA']
sample_names=sort(as.character(sample_names))


# write table:
write.table(sample_names,file=output_file,quote=F,row.names=F,col.names=F)



