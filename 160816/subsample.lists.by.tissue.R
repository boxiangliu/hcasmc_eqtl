#!/usr/bin/Rscript
# 6/1/16
# Joe Davis
# R script to choose which individuals to subsample for in each tissue at each sample size

# Define subsample sizes
subsample.sizes = 52

# Read in list of tissues
tissues = read.table('/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/gtex.v6p.eqtl.tissues.txt', sep = '\t', header = F, stringsAsFactors = F)[, 1]


# subsample:
set.seed(1)
for(i in 1:length(tissues)){
	tissue.samples = t(read.table(paste('/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/v6p/v6p_fastQTL_FOR_QC_ONLY/', tissues[i], 
		'_Analysis.v6p.FOR_QC_ONLY.normalized.expr.bed', sep = ''), nrows = 1, header = F, stringsAsFactors = F, comment.char = '?'))[-c(1:4)]
	subsample = sort(sample(tissue.samples, subsample.sizes))
	out_dir=paste0('../processed_data/160816/', 'subsampling/', tissues[i],'/')
	if (!dir.exists(out_dir)) dir.create(out_dir,recursive=T)
	out_file=paste(out_dir, tissues[i], '_subsample.', subsample.sizes, '.txt', sep = '')
	message(out_file)
	write.table(subsample,out_file, col.names = F, row.names = F, quote = F)
}

