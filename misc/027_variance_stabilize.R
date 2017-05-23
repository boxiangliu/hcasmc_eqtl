#!/usr/bin/env Rscript 
# bosh liu
# 2016/04/27
# durga
# variance stabilize the expression matrix

# read expression matrix: 
count_matrix_filename = '../processed_data/count_matrix.tsv'
count_matrix = fread(count_matrix_filename)


# convert count_matrix from data.table to matrix: 
temp = as.data.frame(count_matrix);
rownames(temp) = temp$gene_id;
temp$gene_id = NULL;
temp = as.matrix(temp);
count_matrix = temp


# remove genes with zero total counts: 
temp = count_matrix[rowSums(count_matrix)!=0,]
count_matrix = temp


# load DESeq2: 
library(DESeq2)


# create DESeq DDS object from count matrix: 
colData = data.frame(intercept = rep(1, ncol(count_matrix)))
colData = rep('HiSeq,length250',ncol(count_matrix))
for (i in 1:ncol(count_matrix)){
	name = colnames(count_matrix)[i]
	new_colData=colData[i]
	if (grepl('dase',name)) new_colData = 'HiSeq,length200'
	if (grepl('CA2305',name)) new_colData = 'NextSeq,length150'
	if (grepl('90702_Nextseq',name)) new_colData = 'NextSeq,length300'
	colData[i]=new_colData
}
temp = data.frame(format=colData)
colData = temp
colnames(dds) = colnames(count_matrix)
dds = DESeqDataSetFromMatrix(countData = count_matrix, colData=colData, design=~format)
dds

# log transformation: 
nt = normTransform(dds) # log2 transformation after library size correction.

# variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds)


# rlog transformation: 
rld <- rlog(dds) 


# load library vsn: 
library(vsn)


# make figure output directory: 
dir.create('../figures/027_variance_stabilize')


# make mean-sd plot for raw counts:
pdf('../figures/027_variance_stabilize/mean_sd.raw_counts.pdf')
meanSdPlot(assay(dds))
dev.off()

# make mean-sd plot for variance stabilized counts: 
pdf('../figures/027_variance_stabilize/mean_sd.variance_stabilized.pdf')
meanSdPlot(assay(vsd))
dev.off()


# make mean-sd plot for rlog transformed counts:
pdf('../figures/027_variance_stabilize/mean_sd.rlog.pdf')
meanSdPlot(assay(rld))
dev.off()


# keep genes with > 20 counts in more than 50% of samples:
for (cutoff in c(20,50)){
	temp = dds[rowSums(counts(dds) > cutoff) > dim(dds)[2]/2,]
	dds2 = temp


	# log transformation: 
	nt2 = normTransform(dds2) # log2 transformation after library size correction.


	# variance stabilizing transformation
	vsd2 <- varianceStabilizingTransformation(dds2)


	# rlog transformation: 
	rld2 <- rlog(dds2) 


	# make mean-sd plot for raw counts:
	pdf(sprintf('../figures/027_variance_stabilize/mean_sd.filtered_%s.raw_counts.pdf',cutoff))
	meanSdPlot(assay(dds2))
	dev.off()

	# make mean-sd plot for variance stabilized counts: 
	pdf(sprintf('../figures/027_variance_stabilize/mean_sd.filtered_%s.variance_stabilized.pdf',cutoff))
	meanSdPlot(assay(vsd2))
	dev.off()


	# make mean-sd plot for rlog transformed counts:
	pdf(sprintf('../figures/027_variance_stabilize/mean_sd.filtered_%s.rlog.pdf',cutoff))
	meanSdPlot(assay(rld2))
	dev.off()
} 


# calculate PCs using variance stabilized data:
dds_filtered = dds[rowSums(counts(dds) > 50) > dim(dds)[2]/2,]
vsd_filtered = varianceStabilizingTransformation(dds_filtered)


# make PCA plot using supplied function:
pdf('../figures/027_variance_stabilize/PCA.pdf')
plotPCA(vsd_filtered,intgroup='format')
dev.off()

save.image('../r_environments/027_variance_stabilize.RData')
# make PCA plot using custom function:
# pcs = prcomp(t(assay(vsd_filtered)),center=T,scale=T)


# # make scree plot:
# pdf('../figures/027_variance_stabilize/scree_plot.pdf')
# plot(pcs,type='l', main = 'scree plot')
# dev.off()


# # calculate percent variance explained by each PC:
# pct_variance = (pcs$sdev)^2/sum((pcs$sdev)^2)


# # plot the first 2 PCs:
# plot(pcs$x[,1:2],xlab=sprintf('PC1 (%0.1f%%)',100*pct_variance[1]), ylab=sprintf('PC2 (%0.1f%%)',100*pct_variance[2]))
# text(pcs$x[,1:2], sample_short)
# names(assay(vsd_filtered))
# head(assay(dds))


