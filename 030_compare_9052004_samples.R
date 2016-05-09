#!/usr/bin/env Rscript
# 2016/05/07
# bosh liu
# durga
# compare 9052004 and 9052004_dase


# libraries: 
library(GGally)
source('utils.R')

# paths: 
count_matrix_filename = '../processed_data/030_read_htseq_count/count_matrix.tsv'
figure_path='../figures/030_compare_9052004_samples'
if (!dir.exists(figure_path)){dir.create(figure_path)}

# read expression matrix: 
count_matrix = fread(count_matrix_filename,header=T)


# convert count_matrix from data.table to matrix: 
temp = as.data.frame(count_matrix);
rownames(temp) = temp$gene_id;
temp$gene_id = NULL;
temp = as.matrix(temp);
count_matrix = temp


# remove genes with zero total counts: 
temp = count_matrix[rowSums(count_matrix)!=0,]
count_matrix = temp


# remove special features:
special_features=c('__alignment_not_unique','__ambiguous','__no_feature')
count_matrix=count_matrix[!rownames(count_matrix)%in%special_features,]


# plot 9052004 vs 9052004_dase
# I also included 1020301 as a control sample
pdf(paste(figure_path,'pairs.pdf',sep='/'))
pairs(count_matrix[,c('9052004_dase','9052004','1020301')],log='xy',upper.panel=panel.cor)
dev.off()


