#!/usr/bin/env Rscript 
# bosh liu
# 2016/04/27
# durga
# variance stabilize the expression matrix


# paths: 
count_matrix_filename = '../processed_data/030_read_htseq_count/count_matrix.tsv'


# read expression matrix: 
count_matrix = fread(count_matrix_filename,header=T)


# convert count_matrix from data.table to matrix: 
temp = as.data.frame(count_matrix);
rownames(temp) = temp$gene_id;
temp$gene_id = NULL;
count_matrix = as.matrix(temp);



# remove genes with zero total counts: 
count_matrix = count_matrix[rowSums(count_matrix)!=0,]



# remove special features:
special_features=c('__alignment_not_unique','__ambiguous','__no_feature')
count_matrix=count_matrix[!rownames(count_matrix)%in%special_features,]


# merge counts from 9052004_dase and 9052004: 
colnames(count_matrix)
temp=count_matrix[,c('9052004','9052004_dase')]
temp=rowSums(temp)
count_matrix[,'9052004']=temp
count_matrix=count_matrix[,colnames(count_matrix)!='9052004_dase']


# load DESeq2: 
library(DESeq2)


# create DESeq DDS object from count matrix: 
colData = data.frame(intercept = rep(1, ncol(count_matrix)))
colData = rep('HiSeq,length250',ncol(count_matrix))
for (i in 1:ncol(count_matrix)){
	name = colnames(count_matrix)[i]
	new_colData=colData[i]
	if (grepl('dase',name)) new_colData = 'HiSeq,length200'
	if (grepl('2305',name)) new_colData = 'NextSeq,length150'
	if (grepl('90702_Nextseq',name)) new_colData = 'NextSeq,length300'
	colData[i]=new_colData
}
colData = data.frame(format=colData)
dds = DESeqDataSetFromMatrix(countData = count_matrix, colData=colData, design=~format)
colnames(dds) = colnames(count_matrix)


# estimate size factors:
dds=estimateSizeFactors(dds,'ratio')
counts_size_corrected=counts(dds,normalized=T)


# save size factor corrected counts:
saveRDS(counts_size_corrected,file='../processed_data/030_variance_stabilize/counts_size_corrected.rds')


# variance stabilizing transformation on all genes:
vsd = varianceStabilizingTransformation(dds)


# keep genes with more than 50 reads in more than half of the samples:
dds_filtered = dds[rowSums(counts(dds) > 50) > dim(dds)[2]/2,]


# variance stabilize on filtered genes:
vsd_filtered = varianceStabilizingTransformation(dds_filtered)


# load library vsn: 
library(vsn)


# make figure output directory: 
dir.create('../figures/030_variance_stabilize')


# make mean-sd plot for variance stabilized counts: 
pdf('../figures/030_variance_stabilize/mean_sd.variance_stabilized.pdf')
meanSdPlot(assay(vsd_filtered))
dev.off()


# save variance stablized counts: 
if (!dir.exists('../processed_data/030_variance_stabilize/')) {dir.create('../processed_data/030_variance_stabilize/')}
saveRDS(vsd,file='../processed_data/030_variance_stabilize/vsd.rds')
saveRDS(vsd_filtered,file='../processed_data/030_variance_stabilize/vsd_filtered.rds')


# make PCA plot using supplied function:
pdf('../figures/030_variance_stabilize/PCA.pdf')
plotPCA(vsd_filtered,intgroup='format')
dev.off()




