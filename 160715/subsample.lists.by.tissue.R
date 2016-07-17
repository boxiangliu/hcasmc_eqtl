#!/usr/bin/Rscript
# 6/1/16
# Joe Davis
# R script to choose which individuals to subsample for in each tissue at each sample size

# Load required packages
require(ggplot2)
require(reshape2)
require(plyr)

#------------ FUNCTIONS

#------------ MAIN

# Define subsample sizes
subsample.sizes = rev(c(70, seq(100, 350, 25)))

# Read in list of tissues
tissues = read.table('/users/joed3/GTExCisEqtls/data/gtex.v6p.eqtl.tissues.txt', sep = '\t', header = F, stringsAsFactors = F)[, 1]

# For each tissue, read in the list of individuals and subsample that list down until the smallest subsample size is reached
# Output each list to the correct tissue directory for use in the FastQTL pipeline
for(i in 1:length(tissues)){
	tissue.samples = t(read.table(paste('/users/joed3/GTExCisEqtls/data/V6P/v6p_fastQTL_FOR_QC_ONLY/', tissues[i], 
		'_Analysis.v6p.FOR_QC_ONLY.normalized.expr.bed', sep = ''), nrows = 1, header = F, stringsAsFactors = F, comment.char = '?'))[-c(1:4)]
	subsample = tissue.samples
	write.table(subsample, paste('/users/joed3/GTExCisEqtls/data/', 'subsampling/', tissues[i], '/', tissues[i], '_subsample.full.txt', sep = ''), sep = '\t', col.names = F, row.names = F, quote = F)
	for(j in 1:length(subsample.sizes)){
		if(length(subsample) >= subsample.sizes[j]){
			subsample = sort(sample(subsample, subsample.sizes[j]))
			write.table(subsample, paste('/users/joed3/GTExCisEqtls/data/', 'subsampling/', tissues[i], '/', tissues[i], '_subsample.', subsample.sizes[j], '.txt', sep = ''), col.names = F, row.names = F, quote = F)
		}
	}
}

