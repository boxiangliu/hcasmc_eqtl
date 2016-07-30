# run ruvseq


# install library:
# source("https://bioconductor.org/biocLite.R")
# biocLite('Biostrings')
# biocLite('Rsamtools')
# biocLite('GenomicAlignments')
# biocLite('ShortRead')
# biocLite('EDASeq')
# biocLite('RUVSeq')
library(RUVSeq)
library(DESeq2)
library("BiocParallel")
register(MulticoreParam(10))


# function: 
#' get OLS residual of counts onto covs:
#' @param counts (matrix) response variables
#' @param covs (matrix) covariants
getResiduals=function(counts,covs){
	residual_mat=matrix(0,nrow=nrow(counts),ncol=ncol(counts),dimnames=list(rownames(counts),colnames(counts)))
	stopifnot(colnames(counts)==rownames(covs))
	M=nrow(counts)
	for (i in 1:M){
		fit=lm(counts[i,]~as.matrix(covs))
		residuals=fit$residuals
		residual_mat[i,]=residuals
	}
	return(residual_mat)
}


# load DESeq2 dataset: 
load('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/find_housekeeping_genes.RData')


# construct seq expression set:
counts=counts(dds)
colData=colData(dds)
tissue=colData$tissue
colnames(counts)=colData$sample
set=newSeqExpressionSet(counts,phenoData=data.frame(tissue=tissue,row.names=colData$sample))


# make RLE and PCA plots: 
library(RColorBrewer)
qualitative_color=brewer.pal(12, "Paired")
pal=colorRampPalette(qualitative_color)
pal=pal(length(unique(colData$tissue)))
color=as.numeric(tissue)

pdf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/rle_and_pca_before_ruv.pdf')
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=pal[color])
plotPCA(set, col=pal[tissue], cex=1.2)
dev.off()


# run RUVg:
hk_genes=read.table('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/hk_genes.txt',header=T)

for (k in seq(5,50,5)) {
	set1=RUVg(set, as.character(hk_genes$Name), k=k)
	pdf(sprintf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/rle_and_pca_after_ruv.%s.pdf',k))
	plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=pal[color])
	plotPCA(set1, col=pal[color], cex=1.2,labels=F)
	dev.off()
}

# Determine if the unwanted factors from a smaller k is the same as 
# the first k factors of a larger K. 
# And yes they are the same. 
test10=RUVg(set, as.character(hk_genes$Name), k=10)
test20=RUVg(set, as.character(hk_genes$Name), k=20)
stopifnot(test10@phenoData@data$W_1==test20@phenoData@data$W_1)
stopifnot(test10@phenoData@data$W_10==test20@phenoData@data$W_10)


# obtain residuals: 
set1=RUVg(set, as.character(hk_genes$Name), k=50)
counts=counts(set1)
covs=pData(set1)
residuals=getResiduals(counts,covs[,2:21]) # using the first 20 factors


# output residuals: 
write.table(residuals,file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/residuals.txt',quote=F,sep='\t',row.names=T,col.names=T)