#!/usr/bin/env Rscript
# bosh liu
# 2016/05/15
# durga
# make sample file for individuals in each ethnicity 


# library:
library(XLConnect)


# paths:
input_file='../processed_data/rna_wgs_match.reduced_050616.xlsx'
output_dir='../processed_data/160526/detect_sample_contamination_model_based'

# read input:
input=readWorksheet(loadWorkbook(input_file),sheet=4)


for (ethnicity in c('AA','Hispanic','Caucasian','Asian')){
	# select ethnic samples:
	sample_names=input[input$Genomic_Ethnicity==ethnicity,'DNA']
	sample_names=sort(as.character(sample_names))

	# write sample names to table:
	output_file=paste0(output_dir,"/",ethnicity,".txt")
	write.table(sample_names,file=output_file,quote=F,row.names=F,col.names=F)
}


