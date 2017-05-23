#!/usr/bin/env Rscript
# bosh liu
# 2016/04/21
# durga
# variance stabilize the count matrix
library(cowplot)
library(tools)


count_dir='../data/rnaseq/expression/working_set/'
sample_files=list.files(count_dir,pattern='Aligned.out.sorted.gene.count',recursive=T,full.name=T)
table_path = '../processed_data/030_read_htseq_count/count_matrix.tsv'
if (!dir.exists(dirname(table_path))) {dir.create(dirname(table_path),recursive=T)}

# generate count matrix
count_matrix = data.table()
for (sample_file in sample_files){
	message(sample_file)
	sample_name = basename(dirname(sample_file))
	gene_count = fread(sample_file)
	setnames(gene_count,c('gene_id','count'))
	gene_count[,sample := sample_name]
	count_matrix = rbind(count_matrix, gene_count)
}

# save count matrix:
temp = data.table::dcast(count_matrix, gene_id ~ sample, value.var = 'count')
count_matrix = temp
write.table(count_matrix, file = table_path, quote=F, sep = '\t', row.names=F, col.names=T)