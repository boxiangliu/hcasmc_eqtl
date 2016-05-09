#!/usr/bin/env Rscript 
# bosh liu
# 2016/04/26
# durga
# perform PCA analysis on expression matrix

# libraries:
library(DESeq2)
library(XLConnect)
source('utils.R')

# paths: 
figure_path='../figures/030_expression_PCA/'
if (!dir.exists(figure_path)) {dir.create(figure_path)}

# read variance stabilized counts: 
vsd=readRDS('../processed_data/030_variance_stabilize/vsd.rds')
vsd_filtered=readRDS('../processed_data/030_variance_stabilize/vsd_filtered.rds')


# perform PCA on vsd_filtered:
pcs = prcomp(t(assay(vsd_filtered)),center=T,scale=T)


# make scree plot:
pdf(paste0(figure_path,'scree_plot.pdf'))
plot((pcs$sdev)^2,main='Scree plot',ylab='Variance',type='b')
dev.off()


# calculate percent variance explained by each PC:
pct_variance = (pcs$sdev)^2/sum((pcs$sdev)^2)


# plot PC1 and PC2:
pdf(paste0(figure_path,'pca.pdf'))
plot(pcs$x[,1:2],xlab=sprintf('PC1 (%0.1f%%)',100*pct_variance[1]), ylab=sprintf('PC2 (%0.1f%%)',100*pct_variance[2]),xlim=c(-200,100))
text(pcs$x[,1:2],rownames(pcs$x))


# plot PC1 and PC3:
plot(pcs$x[,c(1,3)],xlab=sprintf('PC1 (%0.1f%%)',100*pct_variance[1]), ylab=sprintf('PC3 (%0.1f%%)',100*pct_variance[3]),xlim=c(-200,100))
text(pcs$x[,c(1,3)],rownames(pcs$x))


# plot PC2 and PC3:
plot(pcs$x[,c(2,3)],xlab=sprintf('PC2 (%0.1f%%)',100*pct_variance[2]), ylab=sprintf('PC3 (%0.1f%%)',100*pct_variance[3]),xlim=c(-200,100))
text(pcs$x[,c(2,3)],rownames(pcs$x))
dev.off()



# read in covariates table: 
sample_sheet=readWorksheet(loadWorkbook("../processed_data/rna_wgs_match.reduced_050616.xlsx"),sheet=2)
sample_sheet=as.data.table(sample_sheet)


# remove 9072004_dase:
sample_sheet=sample_sheet[RNA!='9052004_dase']


# cast DNA and RNA columns of sample_sheet into character: 
sample_sheet$DNA=as.character(sample_sheet$DNA)
sample_sheet$RNA=as.character(sample_sheet$RNA)


# cast age into integers: 
sample_sheet$Age=as.numeric(sample_sheet$Age)


# sort sample_sheet according to RNA column: 
sample_sheet=sample_sheet[order(sample_sheet$RNA),]



# construct design matrix for Sex: 
design=model.matrix(~Sex+Manufacturer+Age_Imputed+AFR+AMR+EAS+EUR,sample_sheet)
design=design[,-1] # remove intercept term


# correlation between PCs and covariates:
pc_cor=t(cor(pcs$x[,1:10],design))
pc_cor=as.data.table(pc_cor,keep.rownames=T)
setnames(pc_cor,'rn','covariates')


# write correlation to table:
table_dir='../processed_data/030_expression_PCA'
if(!dir.exists(table_dir)){dir.create(table_dir)}
write.table(pc_cor, file=paste(table_dir,'cor_top_10_pcs.tsv',sep="/"),col.names=T,row.names=F,quote=F,sep='\t')


# I manually highlighted the most correlated covariate: 
abs_pc_cor=abs(pc_cor)
max.col(t(abs_pc_cor))


