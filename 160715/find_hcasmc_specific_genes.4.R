#!/usr/bin/env Rscript
# boxiang liu
# durga
# differential expression 

# library:
source('/srv/persistent/bliu2/HCASMC_eQTL/scripts/utils.R')
library("BiocParallel")
register(MulticoreParam(10))
library('DESeq2')
library('gap')
library('gplots')

# functions:
select_top_genes=function(norm_counts,resOrdered,col_data,tissues,n_genes=100,up=TRUE){
	# order the rows of count matrix:
	if (up==TRUE){
		idx=match(resOrdered[log2FoldChange>=0,gene_id],rownames(norm_counts))
	} else {
		idx=match(resOrdered[log2FoldChange<=0,gene_id],rownames(norm_counts))
	}
	norm_counts_ordered=norm_counts[idx,]


	# extract top DE genes:
	gene_select=1:n_genes
	tissue_select=which(col_data$tissue%in%c('HCASMC',tissues))
	top_genes=norm_counts_ordered[gene_select,tissue_select]


	# give sensible names to top_genes: 
	rownames(top_genes)=resOrdered[gene_select,gene_name]
	colnames(top_genes)=col_data$tissue[tissue_select]


	# create color for column color strip: 
	tissue_colors=ifelse(col_data$tissue[tissue_select]=='HCASMC',"#FF0000","#0000FF")


	# return: 
	return(list(top_genes=top_genes,tissue_colors=tissue_colors))
}


# read input:
hcasmc_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/combined.count'
hcasmc=read.table(hcasmc_file,header=T,check.names=F)
hcasmc_sub=subsample_gct(hcasmc,10)
res=decompose_gct(hcasmc_sub)


# extract gene_id and gene_name from hcascm: 
row_data=hcasmc_sub%>%select(Name,Description)
setnames(row_data,c('gene_id','gene_name'))


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


# create DESeq dataset: 
dds=DESeqDataSetFromMatrix(countData = as.matrix(master_count),colData=master_col_data,design=~tissue+0)
dds=dds[rowSums(counts(dds))>=10,]



# run DESeq:
dds=DESeq(dds,betaPrior=FALSE,parallel=T)


#### HCASMC vs all GTEx tissue:
# get contrast result: 
res=results(dds, contrast = c(50,rep(-1,50))/51)


# plot the distribution of p-values: 
pdf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/DE_p_value_distribution.3.pdf')
hist(res$pvalue,main='DE p-value',xlab='p-value')
qqunif(res$pvalue)
dev.off()


# cast result into data.table:
res2=data.table(as.data.frame(res),keep.rownames=T)
setnames(res2,'rn','gene_id')


# add gene_name column to result: 
res3=merge(res2,row_data,by='gene_id',all.x=T)
setcolorder(res3,c(1,ncol(res3),2:(ncol(res3)-1)))


# sort result by pvalue:
resOrdered=res3%>%arrange(pvalue)


# look at the top few genes: 
resOrdered[log2FoldChange>0,.(gene_name,gene_id)][1:5]
#    gene_name            gene_id
# 1:     RPL21  ENSG00000122026.6
# 2:     RPS26  ENSG00000197728.5
# 3:    RPL36A  ENSG00000241343.5
# 4:     YIPF5 ENSG00000145817.12
# 5:    DCBLD2 ENSG00000057019.11
# the first 3 genes are ribosomal genes; YIPF5 transports between ER and Golgi;
# DCBLD2, also named Endothelial and smooth muscle-derived neuropilin-like protein, regulates platelet-derived growth factor
# signaling in human vascular smooth muscle cells by modulating receptor ubiquitination.

# extract normalized counts: 
norm_counts=counts(dds,normalized=T)


# extract column data:
col_data=colData(dds)


# extract the gene_name and gene_id of the top 5 genes: 
top5up=resOrdered[log2FoldChange>0,.(gene_name,gene_id)][1:5]
top5dn=resOrdered[log2FoldChange<0,.(gene_name,gene_id)][1:5]


# plot the top 5 upregulated genes:
pdf('../figures/160715/DE_up_examples.pdf',onefile=T)
for (i in 1:nrow(top5up)){
	gene_id=top5up$gene_id[i]
	gene_name=top5up$gene_name[i]
	gene_counts=norm_counts[rownames(norm_counts)==gene_id]
	to_plot=data.frame(tissue=col_data$tissue,counts=gene_counts)
	p=ggplot(to_plot,aes(reorder(tissue, -counts, FUN=median),counts))+geom_boxplot()+scale_y_log10()+theme_bw()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+ylab('Normalized counts')+xlab('Tissue')+ggtitle(gene_name)
	print(p)
}
dev.off()


# plot the top 5 down-regulated genes:
pdf('../figures/160715/DE_down_examples.pdf',onefile=T)
for (i in 1:nrow(top5dn)){
	gene_id=top5dn$gene_id[i]
	gene_name=top5dn$gene_name[i]
	gene_counts=norm_counts[rownames(norm_counts)==gene_id]
	to_plot=data.frame(tissue=col_data$tissue,counts=gene_counts)
	p=ggplot(to_plot,aes(reorder(tissue, -counts, FUN=median),counts))+geom_boxplot()+scale_y_log10()+theme_bw()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+ylab('Normalized counts')+xlab('Tissue')+ggtitle(gene_name)
	print(p)
}
dev.off()



# order the rows of count matrix: 
idx=match(resOrdered$gene_id,rownames(norm_counts))
norm_counts_ordered=norm_counts[idx,]


# extract top 5000 DE genes:
gene_select=1:200
tissue_select=c(1:10,seq(11,510,10))
top_genes=norm_counts_ordered[gene_select,tissue_select]


# give sensible names to top_genes: 
rownames(top_genes)=resOrdered[gene_select,gene_name]
colnames(top_genes)=col_data$tissue[tissue_select]


# create color for column color strip: 
tissue_colors=ifelse(col_data$tissue[tissue_select]=='HCASMC',"#FF0000","#0000FF")
resOrdered[log2FoldChange>0,][1:10]


# make heatmap of DE genes:
dist2=function(x) {pearson=cor(t(x));distance=as.dist(1-pearson);return(distance)}
pdf('../figures/160715/heatmap.contrast_model.pdf')
heatmap.2(top_genes, distfun=dist2,col=topo.colors(75), scale="row", ColSideColors=tissue_colors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
dev.off()




# output rank:
output=resOrdered[,.(gene_name,stat)]
write.table(output,file='../processed_data/160715/hcasmc_vs_gtex_all.rnk',quote=F,row.names=F,col.names=F,sep='\t')


#### HCASMC vs GTEx artery tissue:
# setup: 
artery_tissues=c('Artery_Aorta','Artery_Coronary','Artery_Tibial')
contrast=c(length(artery_tissues),rep(-1,50))/(length(artery_tissues)+1)


# perform contrast testing procedure: 
res_artery=results(dds,contrast=ifelse(unique(col_data$tissue)%in%c('HCASMC',artery_tissues),contrast,0))


# reorder rows by p-value: 
resOrdered_artery=format_result(res_artery,row_data)


# save: 
save(res_artery,resOrdered_artery,file='../processed_data/160715/hcasmc_vs_gtex_artery.RData')


# select top upregulated genes: 
temp=select_top_genes(norm_counts,resOrdered_artery,col_data,tissues=artery_tissues,n_genes=100,up=TRUE)
top_genes_up=temp[[1]]
tissue_colors=temp[[2]]


# make heatmap of upregulated genes:
pdf('../figures/160715/heatmap.contrast_model.hcasmc_vs_GTEx_artery.upregulation.pdf')
heatmap.2(top_genes_up,Colv=F,Rowv=F,dendrogram='none',col=topo.colors(75), scale="row", ColSideColors=tissue_colors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
dev.off()


# select top downregulated genes:  
temp=select_top_genes(norm_counts,resOrdered_artery,col_data,tissues=artery_tissues,n_genes=100,up=FALSE)
top_genes_down=temp[[1]]
tissue_colors=temp[[2]]


# make heatmap of down regulated genes: 
pdf('../figures/160715/heatmap.contrast_model.hcasmc_vs_GTEx_artery.downregulation.pdf')
heatmap.2(top_genes_down,Colv=F,Rowv=F,dendrogram='none',col=topo.colors(75), scale="row", ColSideColors=tissue_colors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
dev.off()


# output rank: 
output=resOrdered_artery[,.(gene_name,stat)]
write.table(output,file='../processed_data/160715/hcasmc_vs_gtex_artery.rnk',quote=F,row.names=F,col.names=F,sep='\t')


#### smooth muscle: 
# setup:
target_tissue=c('Bladder','Colon_Sigmoid','Colon_Transverse','Esophagus_Gastroesophageal_Junction','Esophagus_Mucosa', 'Esophagus_Muscularis','Small_Intestine_Terminal_Ileum','Stomach','Uterus')
contrast=c(length(target_tissue),rep(-1,50))/(length(target_tissue)+1)


# perform contrast testing procedure: 
res=results(dds,contrast=ifelse(unique(col_data$tissue)%in%c('HCASMC',target_tissue),contrast,0))


# reorder rows by p-value: 
resOrdered=format_result(res,row_data)


# save: 
save(res,resOrdered,file='../processed_data/160715/hcasmc_vs_gtex_smooth_muscle.RData')


# select top upregulated genes: 
temp=select_top_genes(norm_counts,resOrdered,col_data,tissues=target_tissue,n_genes=100,up=TRUE)
top_genes_up=temp[[1]]
tissue_colors=temp[[2]]


# make heatmap of upregulated genes:
pdf('../figures/160715/heatmap.contrast_model.hcasmc_vs_GTEx_smooth_muscle.upregulation.pdf')
heatmap.2(top_genes_up,Colv=F,Rowv=F,dendrogram='none',col=topo.colors(75), scale="row", ColSideColors=tissue_colors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
dev.off()


# select top downregulated genes:  
temp=select_top_genes(norm_counts,resOrdered,col_data,tissues=target_tissue,n_genes=100,up=FALSE)
top_genes_down=temp[[1]]
tissue_colors=temp[[2]]


# make heatmap of down regulated genes: 
pdf('../figures/160715/heatmap.contrast_model.hcasmc_vs_GTEx_smooth_muscle.downregulation.pdf')
heatmap.2(top_genes_down,Colv=F,Rowv=F,dendrogram='none',col=topo.colors(75), scale="row", ColSideColors=tissue_colors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
dev.off()


# output rank: 
output=resOrdered[,.(gene_name,stat)]
write.table(output,file='../processed_data/160715/hcasmc_vs_gtex_smooth_muscle.rnk',quote=F,row.names=F,col.names=F,sep='\t')



#### heart: 
# setup:
target_tissue=c('Heart_Atrial_Appendage','Heart_Left_Ventricle')
contrast=c(length(target_tissue),rep(-1,50))/(length(target_tissue)+1)


# perform contrast testing procedure: 
res=results(dds,contrast=ifelse(unique(col_data$tissue)%in%c('HCASMC',target_tissue),contrast,0))


# reorder rows by p-value: 
resOrdered=format_result(res,row_data)


# save: 
save(res,resOrdered,file='../processed_data/160715/hcasmc_vs_gtex_heart.RData')


# select top upregulated genes: 
temp=select_top_genes(norm_counts,resOrdered,col_data,tissues=target_tissue,n_genes=100,up=TRUE)
top_genes_up=temp[[1]]
tissue_colors=temp[[2]]


# make heatmap of upregulated genes:
pdf('../figures/160715/heatmap.contrast_model.hcasmc_vs_GTEx_heart.upregulation.pdf')
heatmap.2(top_genes_up,Colv=F,Rowv=F,dendrogram='none',col=topo.colors(75), scale="row", ColSideColors=tissue_colors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
dev.off()


# select top downregulated genes:  
temp=select_top_genes(norm_counts,resOrdered,col_data,tissues=target_tissue,n_genes=100,up=FALSE)
top_genes_down=temp[[1]]
tissue_colors=temp[[2]]


# make heatmap of down regulated genes: 
pdf('../figures/160715/heatmap.contrast_model.hcasmc_vs_GTEx_heart.downregulation.pdf')
heatmap.2(top_genes_down,Colv=F,Rowv=F,dendrogram='none',col=topo.colors(75), scale="row", ColSideColors=tissue_colors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
dev.off()


# output rank: 
output=resOrdered[,.(gene_name,stat)]
write.table(output,file='../processed_data/160715/hcasmc_vs_gtex_heart.rnk',quote=F,row.names=F,col.names=F,sep='\t')



#### heart: 
# setup:
target_tissue=c('Cells_Transformed_fibroblasts')
contrast=c(length(target_tissue),rep(-1,50))/(length(target_tissue)+1)


# perform contrast testing procedure: 
res=results(dds,contrast=ifelse(unique(col_data$tissue)%in%c('HCASMC',target_tissue),contrast,0))


# reorder rows by p-value: 
resOrdered=format_result(res,row_data)


# save: 
save(res,resOrdered,file='../processed_data/160715/hcasmc_vs_gtex_fibroblast.RData')


# select top upregulated genes: 
temp=select_top_genes(norm_counts,resOrdered,col_data,tissues=target_tissue,n_genes=100,up=TRUE)
top_genes_up=temp[[1]]
tissue_colors=temp[[2]]


# make heatmap of upregulated genes:
pdf('../figures/160715/heatmap.contrast_model.hcasmc_vs_GTEx_fibroblast.upregulation.pdf')
heatmap.2(top_genes_up,Colv=F,Rowv=F,dendrogram='none',col=topo.colors(75), scale="row", ColSideColors=tissue_colors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
dev.off()


# select top downregulated genes:  
temp=select_top_genes(norm_counts,resOrdered,col_data,tissues=target_tissue,n_genes=100,up=FALSE)
top_genes_down=temp[[1]]
tissue_colors=temp[[2]]


# make heatmap of down regulated genes: 
pdf('../figures/160715/heatmap.contrast_model.hcasmc_vs_GTEx_fibroblast.downregulation.pdf')
heatmap.2(top_genes_down,Colv=F,Rowv=F,dendrogram='none',col=topo.colors(75), scale="row", ColSideColors=tissue_colors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
dev.off()


# output rank: 
output=resOrdered[,.(gene_name,stat)]
write.table(output,file='../processed_data/160715/hcasmc_vs_gtex_fibroblast.rnk',quote=F,row.names=F,col.names=F,sep='\t')