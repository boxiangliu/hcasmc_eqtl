#!/usr/bin/env Rscript 
# bosh liu
# 2016/04/20
# durga
# concatenate output from impute2

args = commandArgs(T)
# wd = '../data/joint/imputation/'
wd = args[1]
for (chrom in seq(1,22)){
	pattern = sprintf('recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.chr%s.phased.imputed.[0-9]+_[0-9]+$',chrom)
	files = list.files(path = wd, pattern = pattern, full.names = T)
	region = str_split_fixed(files, pattern = '\\.', n = 15)[,15]
	start = str_split_fixed(region, pattern = '_', n = 2)[,1]
	start = as.integer(start)
	sorted_files = files[order(start)]
	cmd = paste('cat', paste(sorted_files, collapse = ' '), sprintf("> %s/chr%s.impute2", wd, chrom))
	system(cmd) 

	# 
	# pattern = sprintf('recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.chr%s.phased.imputed.[0-9]+_[0-9]+_info$',chrom)
	# files = list.files(path = wd, pattern = pattern, full.names = T)
	# region = str_split_fixed(files, pattern = '\\.', n = 15)[,15]
	# start = str_split_fixed(region, pattern = '_', n = 2)[,1]
	# start = as.integer(start)
	# sorted_files = files[order(start)]
	# cmd = paste('cat', paste(sorted_files, collapse = ' '), sprintf("> %s/chr%s.impute2_info", wd, chrom))
	# system(cmd) 
}



