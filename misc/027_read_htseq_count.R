#!/usr/bin/env Rscript
# bosh liu
# 2016/04/21
# durga
# variance stabilize the count matrix
library(cowplot)
library(tools)

#---- generate QC plots -------
args = commandArgs(T)
sample_list = args[1]
# sample_list = 'read_htseq_count.sample_list.txt'
table_path = args[2]
# table_path = '../processed_data/count_matrix.tsv'

for (stranded in c('forward','reverse')){
	count_dir = sprintf('../data/rnaseq/expression/%s/',stranded)
	sample_files = scan(sample_list, what = 'character', sep = "\n")
	count_matrix = data.table()
	for (sample_file in sample_files){
		if (grepl("^#", sample_file)) {
			message('skipping ',sample_file)
			next()
		}
		message(sample_file)
		sample_path = paste(count_dir,sample_file,sep='/')
		sample_name = dirname(sample_path)
		sample_name = basename(sample_name)
		gene_count = fread(sample_path)
		setnames(gene_count,c('gene_id','count'))
		gene_count[,sample := sample_name]
		count_matrix = rbind(count_matrix, gene_count)
	}

	temp = data.table::dcast(count_matrix, gene_id ~ sample, value.var = 'count')
	count_matrix = temp

	pct_special = data.table(pct=numeric(), feature=character(), sample=character())

	for (sample in colnames(count_matrix)[-1]){
		sample_short = paste(str_split(sample,'_')[[1]][1],str_split(sample,'_')[[1]][2], sep='_')
		for (feature in c('__no_feature','__ambiguous','__too_low_aQual','__not_aligned','__alignment_not_unique')){
			temp = data.table(pct = unlist(count_matrix[gene_id == feature,sample, with=F])/sum(count_matrix[,sample,with=F]), feature = feature, sample = sample_short)
			pct_special = rbind(pct_special, temp)
		}
	}
	p = ggplot(pct_special, aes(x=sample,y=pct,color=feature)) + geom_point() + theme(axis.text.x = element_text(angle = 90,vjust=0.5)) + background_grid(major = "x", minor = "none") + geom_text(aes(label=ifelse(pct>0.2,sample,'')),angle = 90, hjust = 1.1) + ylab('Fraction') + xlab('Sample') + ggtitle(toTitleCase(stranded))
	save_plot(sprintf('../figures/gene_count_%s.pdf',stranded), p, base_height = 10)
}

#------ generate count matrix ----------
count_matrix = data.table()
for (sample_file in sample_files){
	if (grepl("^#", sample_file)) {
		message('skipping ',sample_file)
		next()
	}
	message(sample_file)
	sample_name = dirname(sample_file)
	stranded = ifelse(sample_name == '90702_Nextseq', 'forward', 'reverse')
	count_dir = sprintf('../data/rnaseq/expression/%s/',stranded)
	sample_path = paste(count_dir,sample_file,sep='/')
	gene_count = fread(sample_path)
	setnames(gene_count,c('gene_id','count'))
	gene_count[,sample := sample_name]
	count_matrix = rbind(count_matrix, gene_count)
}

temp = data.table::dcast(count_matrix, gene_id ~ sample, value.var = 'count')
count_matrix = temp
write.table(count_matrix, file = table_path, quote=F, sep = '\t', row.names=F, col.names=T)