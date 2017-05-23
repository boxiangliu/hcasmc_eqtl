#!/usr/bin/env Rscript
# bosh liu
# 2016/04/17
# durga
# plot impute2 concordance table

log_filename = "../data/joint/imputation/recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.chr9.phased.imputed.95000001_100000000_summary"
log = fread(log_filename, skip = 153)
log = log[,.(V6, V7,V8)]

setnames(log,c('interval','pct_called','pct_concordance'))
log$interval = str_replace(log$interval, ']','')

plot(log$pct_called, log$pct_concordance, xlab = "% genotypes called", ylab = '% concordance')

for (chr in 1:22){
	pattern = sprintf('recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.chr%s.phased.imputed.*_summary',chr)
	path = '../data/joint/imputation/'
	files = list.files(path,pattern)
	for (file in files){
		skip = 153
		next_loop=FALSE
		log = tryCatch(fread(paste0(path, file), skip = skip), error = function(e) {message(file); next_loop <<- TRUE; message('last line:'); system(paste0('tail -n2 ',paste0(path, file)))})
		if (next_loop) next()
		log = log[,.(V6, V7,V8)]
		setnames(log,c('interval','pct_called','pct_concordance'))
		log$interval = str_replace(log$interval, ']','')
		tryCatch(points(log$pct_called, log$pct_concordance),error = function(e) {plot(log$pct_called, log$pct_concordance, xlab = "% genotypes called", ylab = '% concordance')})
	}
}