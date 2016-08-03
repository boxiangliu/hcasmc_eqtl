#!/usr/bin/env Rscript
# boxiang liu
# durga
# differential expression 

# library:
source('/srv/persistent/bliu2/HCASMC_eQTL/scripts/utils.R')
library("BiocParallel")
register(MulticoreParam(20))
library('DESeq2')

# functions:
# all moved to utils.R 


# read input:
hcasmc_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/combined.count'
hcasmc=read.table(hcasmc_file,header=T,check.names=F)
hcasmc_sub=subsample_gct(hcasmc,10)
res=decompose_gct(hcasmc_sub)


# inititilize: 
master_count=res$count
master_col_data=data.frame(res$col_data,tissue='HCASMC')
master_row_data=res$row_data


gtex_files=list.files('/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/v6p/subsampling',pattern='*.10.count',recursive=T,full.name=T)
for (gtex_file in gtex_files){
	tissue=str_replace(basename(gtex_file),'_subsample.10.count','')
	message(tissue)
	count=read.table(gtex_file,header=T,check.names=F)

	res=decompose_gct(count)
	count=res$count
	col_data=res$col_data
	row_data=res$row_data
	col_data$tissue=tissue


	# sanity check:
	stopifnot(setequal(row_data$Name,master_row_data$Name))
	stopifnot(length(unique(row_data$Name))==length(master_row_data$Name))


	# reorder the rows of count: 
	idx=match(rownames(master_count),rownames(count))
	count=count[idx,]
	row_data=row_data[idx,]
	stopifnot(row_data$Name==master_row_data$Name)
	stopifnot(row.names(count)==row.names(master_count))


	# append to master: 
	master_col_data=rbind(master_col_data,col_data)
	master_count=data.frame(master_count,count,check.names=F)
}
# head(master_count)

# create DESeq dataset: 
dds=DESeqDataSetFromMatrix(countData = as.matrix(master_count),colData=master_col_data,design=~tissue)
dds=dds[rowSums(counts(dds))>=10,]


# run DESeq:
dds=DESeq(dds,test='LRT',reduced=~1,parallel=T)
res=results(dds)
resOrdered=res[order(res$padj),]


# number of genes with FDR > 0.05, 0.1, 0.5, 0.99:
sum(res$padj>0.05) # 10692
sum(res$padj>0.1) # 10496
sum(res$padj>0.5) # 9895
sum(res$padj>0.99) # 9435


# plot the distribution of p-values: 
pdf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/DE_p_value_distribution.pdf')
hist(res$pvalue,breaks=1000,main='DE p-value',xlab='p-value')
plot(res$pvalue,main='DE p-value',ylab='p-value')
dev.off()


# plot the mean read depth vs p-value:
pdf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/mean_read_depth_vs_pvalue.pdf')
with(res,plot(baseMean,pvalue,log='x'))
dev.off()



# make a boxplot of one of the not differentially expressed gene: 
pdf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/read_counts_across_tissue_constant_genes.pdf')
par(mar=c(17.1,4.1,4.1,2.1))
for (gene in row.names(res[which(res$padj>0.1&res$baseMean>1),])){
	par(las=2)
	data=plotCounts(dds,gene=gene,intgroup='tissue',returnData=T)
	with(data,plot(tissue,count,main=gene))

}
dev.off()


# make boxplot for two differentially expressed genes: 
pdf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/read_counts_across_tissue_for_DE_genes.pdf')
par(mar=c(17.1,4.1,4.1,2.1))
for (i in c('ENSG00000139915.14','ENSG00000114861.14')){
	par(las=2)
	gene=rownames(res[i,])
	data=plotCounts(dds,gene=gene,intgroup='tissue',returnData=T)
	with(data,plot(tissue,count,main=gene))
}
dev.off()


# look at several house keeping genes reported by Levenon et al.
hk_genes=c('C1orf43','CHMP2A','EMC7','GPI','PSMB2','PSMB4','RAB7A','REEP5','SNRPD3','VCP','VPS29')
pdf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/read_counts_across_tissue_for_housing_keeping_genes_by_Levanon.pdf')
par(mar=c(17.1,4.1,4.1,2.1))
for (gene in hcasmc[hcasmc$Description%in%hk_genes,'Name']){
	par(las=2)
	gene_name=hcasmc[hcasmc$Name==gene,'Description']
	data=plotCounts(dds,gene=gene,intgroup='tissue',returnData=T)
	with(data,plot(tissue,count,main=gene))
	mtext(gene_name,las=1)
}
dev.off()


# retrieve size factor corrected counts:
counts=counts(dds,normalized=T)


# perform variance stabilizing transformation: 
vsnnbp=varianceStabilizingTransformation(dds,blind=F)
vsnnbl=varianceStabilizingTransformation(dds,fitType='local',blind=F)
vsnnbm=varianceStabilizingTransformation(dds,fitType='mean',blind=F)


# compare log transformation, sqrt, 2/3 power transformation, and vsn:
pdf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/compare_transformations.pdf')
par(mfrow=c(3,2))
hist(counts[2,],breaks=100,main='Raw count data')
hist(counts[2,]^(1/2),breaks=100, main='sqrt(x)')
hist(counts[2,]^(2/3),breaks=100, main='x^(2/3)')
hist(log2(counts[2,]+1),breaks=100, main='log(x+1)')
hist(assay(vsnnbl)[2,],breaks=100,main='parametric VSN')
hist(assay(vsnnbp)[2,],breaks=100,main='local VSN')
dev.off()


# plot the mean variance relationship:
pdf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/mean_vs_variance.pdf')
par(mfrow=c(1,3))
meanSdPlot(counts)
meanSdPlot(log2(counts+1))
meanSdPlot(assay(vsnnbl))
dev.off()


# plot the variance across tissues for one gene: 
pdf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/expression_by_tissue_raw_count_vs_vsn.pdf',width=12,height=6)
par(mfrow=c(2,2))
temp=data.table(tissue=colData(dds)$tissue,count=counts[2,])
temp[,sd:=sd(count),by='tissue']
temp[,mean:=mean(count),by='tissue']
temp[,count_offset:=count-mean]
with(temp,plot(tissue,count_offset,main='raw count'))
fit=lm(sd~mean,temp)
with(temp,plot(mean,sd,main='raw count'))
abline(fit)


count_vsn=assay(vsnnbl)
temp=data.table(tissue=colData(dds)$tissue,count=count_vsn[2,])
temp[,sd:=sd(count),by='tissue']
temp[,mean:=mean(count),by='tissue']
temp[,count_offset:=count-mean]
with(temp,plot(tissue,count_offset,main='vsn'))
fit=lm(sd~mean,temp)
with(temp,plot(mean,sd,main='vsn'))
abline(fit)
dev.off()



# save: 
save(dds,res,resOrdered,count_vsn,file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/find_housekeeping_genes.RData')

# load: 
# load('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/find_housekeeping_genes.RData')


#### RPKM analysis: 
# read input:
hcasmc_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160519_rpkm/combined.rpkm'
hcasmc=read.table(hcasmc_file,header=T,check.names=F)
hcasmc_sub=subsample_gct(hcasmc,10)
res=decompose_gct(hcasmc_sub)


# inititilize: 
master_count=res$count
master_col_data=data.frame(res$col_data,tissue='HCASMC')
master_row_data=res$row_data


gtex_files=list.files('/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/v6p/subsampling',pattern='*.10.rpkm',recursive=T,full.name=T)
for (gtex_file in gtex_files){
	tissue=str_replace(basename(gtex_file),'_subsample.10.rpkm','')
	message(tissue)
	count=read.table(gtex_file,header=T,check.names=F)

	res=decompose_gct(count)
	count=res$count
	col_data=res$col_data
	row_data=res$row_data
	col_data$tissue=tissue


	# sanity check:
	stopifnot(setequal(row_data$Name,master_row_data$Name))
	stopifnot(length(unique(row_data$Name))==length(master_row_data$Name))


	# reorder the rows of count: 
	idx=match(rownames(master_count),rownames(count))
	count=count[idx,]
	row_data=row_data[idx,]
	stopifnot(row_data$Name==master_row_data$Name)
	stopifnot(row.names(count)==row.names(master_count))


	# append to master: 
	master_col_data=rbind(master_col_data,col_data)
	master_count=data.frame(master_count,count,check.names=F)
}


# merge count table and col_data to facilitate melting:
master_count_t=data.frame(t(master_count))
master_count_t$sample=rownames(master_count_t)
merged=merge(master_col_data,master_count_t,by='sample')


# melt merged data to facilitate filtering:
melted=melt(merged,id=c('sample','tissue'),variable.name='gene_id',value.name='rpkm')
melted=melted%>%group_by(gene_id,tissue)%>%mutate(median=median(rpkm))
melted=melted%>%dplyr::select(tissue,gene_id,median)%>%unique()
melted=as.data.table(melted)


# calculate log median rpkm and sd: 
melted[,logmedian:=log2(median+1)]
melted[,sd:=sd(logmedian),by='gene_id']


# apply sd filter, expression filter, and log fold change filter:
filtered=melted[sd<1,]
filtered[,min:=min(median),by='gene_id']
filtered=filtered[min>0,]
filtered[,avg:=mean(logmedian),by='gene_id']
filtered[,maxdiff:=max(logmedian)-avg,by='gene_id']
filtered=filtered[maxdiff<1,]


# get and save housekeeping gene names:
hk_genes=unique(filtered$gene_id)
hk_genes=hcasmc[hcasmc$Name%in%hk_genes,c('Name','Description')]
write.table(hk_genes,'/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/hk_genes.pdf',quote=F,sep='\t',row.names=F)


# plot a housekeeping gene: 
hk_gene_exp=master_count[rownames(master_count)==as.character(hk_genes[1,1]),]
stopifnot(master_col_data$sample==colnames(hk_gene_exp))
pdf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/read_counts_for_hk_gene.txt')
to_plot=data.frame(tissue=master_col_data$tissue,rpkm=unlist(hk_gene_exp))
ggplot(to_plot,aes(tissue,rpkm))+geom_boxplot()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()


# plot a random non housekeeping gene: 
non_hk_gene_exp=master_count[which(!rownames(master_count)%in%hk_genes[,1])[2],]
pdf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/read_counts_for_non_hk_gene.pdf')
to_plot=data.frame(tissue=master_col_data$tissue,rpkm=unlist(non_hk_gene_exp))
ggplot(to_plot,aes(tissue,rpkm))+geom_boxplot()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()