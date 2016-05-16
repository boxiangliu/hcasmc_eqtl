#!/usr/bin/env Rscript
# bosh liu
# 2016/05/15
# durga
# make sample file for 52 samples in the working set


# library:
library(XLConnect)


# paths:
input_file='../processed_data/rna_wgs_match.reduced_050616.xlsx'
output_file='../data/joint2/sample_list.txt'


# read input:
input=readWorksheet(loadWorkbook(input_file),sheet=2)


# sort unique sample names:
sample_names=unique(input$DNA)
sample_names=sort(as.character(sample_names))


# write sample names to table:
write.table(sample_names,file=output_file,quote=F,row.names=F,col.names=F)
