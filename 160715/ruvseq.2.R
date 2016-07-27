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


#### function: 
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


plotMeanAndVar=function(x){
	median=apply(x,1,median)
	x=x-median
	x=as.data.frame(x)
	x$gene_id=rownames(x)
	x=as.data.table(x)
	melted=melt(x,id.vars='gene_id',variable.name='sample',value.name='log_expr')
	melted[,median:=median(log_expr),by='sample']
	melted[,mean:=mean(log_expr),by='sample']
	melted[,variance:=var(log_expr),by='sample']
	temp=unique(melted[,.(sample,mean,median,variance)])
	to_plot=melt(temp,id.vars='sample',variable.name='statistics',value.name='value')
	p=ggplot(to_plot,aes(as.factor(sample),value,color=statistics))+geom_point(alpha=0.5)+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+xlab('sample')
	return(p)
}


plotMeanAndVarByTissue=function(x,hash){
	median=apply(x,1,median)
	x=x-median
	x=as.data.frame(x)
	x$gene_id=rownames(x)
	x=as.data.table(x)
	melted=melt(x,id.vars='gene_id',variable.name='sample',value.name='log_expr')
	melted=merge(melted,hash,by='sample')
	melted[,median:=median(log_expr),by='tissue']
	melted[,mean:=mean(log_expr),by='tissue']
	melted[,variance:=var(log_expr),by='tissue']
	temp=unique(melted[,.(tissue,mean,median,variance)])
	to_plot=melt(temp,id.vars='tissue',variable.name='statistics',value.name='value')
	p=ggplot(to_plot,aes(as.factor(tissue),value,color=statistics))+geom_point()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+xlab('tissue')+background_grid(major='x')
	return(p)
}


#### main
# load DESeq2 dataset: 
load('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/find_housekeeping_genes.RData')


# construct seq expression set:
counts=counts(dds,normalized=T)
colData=colData(dds)
tissue=colData$tissue
colnames(counts)=colData$sample


# load housekeeping genes:
hk_genes=read.table('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/hk_genes.txt',header=T)


# log transform counts:
log_counts=log(counts+1)


# plot mean and variance by sample:
for (k in seq(5,50,5)){
	set1=RUVg(log_counts, as.character(hk_genes$Name), k=k,isLog=T)
	p=plotMeanAndVar(set1$normalizedCounts)
	p=p+theme(axis.text.x=element_blank())
	save_plot(sprintf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/rle_and_pca_after_ruv.2.%s.pdf',k),p)
}

# plot mean and variance by tissue:
hash=as.data.table(colData[,c('sample','tissue')])
for (k in seq(5,50,5)){
	set1=RUVg(log_counts, as.character(hk_genes$Name), k=k,isLog=T)
	p=plotMeanAndVarByTissue(set1$normalizedCounts,hash)
	save_plot(sprintf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/rle_and_pca_after_ruv.2.%s.by_tissue.pdf',k),p,base_width=10,base_height=10)
}


# obtain residuals: 
set1=RUVg(log_counts, as.character(hk_genes$Name), k=20,isLog=TRUE)
residuals=as.data.frame(set1$normalizedCounts)

# output residuals: 
write.table(residuals,file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/residuals.2.20.txt',quote=F,sep='\t',row.names=TRUE,col.names=TRUE)