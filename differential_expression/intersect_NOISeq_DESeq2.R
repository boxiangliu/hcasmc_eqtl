library(data.table)
library(DESeq2)
library(stringr)
library(dplyr)
library(cowplot)
library(foreach)
library(doMC)
registerDoMC(6)

# Variables:
deseq2_fn='../processed_data/differential_expression/DESeq2/dds.rds'
noiseq_fn='../processed_data/differential_expression/noiseq/noiseq.rds'
tissue_list=c("Artery - Aorta","Artery - Coronary","Artery - Tibial", 'Heart - Atrial Appendage', 'Heart - Left Ventricle','Cells - Transformed fibroblasts')
fig_dir='../figures/differential_expression/intersect_NOISeq_DESeq2/'
out_dir='../processed_data/differential_expression/intersect_NOISeq_DESeq2/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir)}
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# Functions:
read_deseq2_result=function(fn){
	deseq2=readRDS(fn)
	return(deseq2)
}


read_noiseq_result=function(fn){
	noiseq=readRDS(fn)
	return(noiseq)
}


select_tissue=function(x,tissue){
	x[[tissue]]
}


merge_noiseq_and_deseq2=function(noiseq,deseq2,tissue){
	message('INFO - getting deseq2 results..')
	deseq2_res=as.data.frame(results(deseq2,contrast=c('tissue','HCASMC',tissue))[,c('log2FoldChange','padj')])
	setnames(deseq2_res,'log2FoldChange','deseq2_log2FC')

	message('INFO - getting noiseq results..')
	noiseq_res=noiseq@results[[1]][,c('prob','log2FC')]
	numerator=str_split_fixed(noiseq@comparison,' - ',2)[,1]
	if (numerator!='HCASMC'){
		noiseq_res$log2FC=-noiseq_res$log2FC
	}
	setnames(noiseq_res,'log2FC','noiseq_log2FC')

	noiseq_res=noiseq_res[match(rownames(deseq2_res),rownames(noiseq_res)),]
	stopifnot(all(rownames(deseq2_res)==rownames(noiseq_res)))

	message('INFO - merging results.. ')
	merged=cbind(deseq2_res,noiseq_res)
	stopifnot(cor(merged$noiseq_log2FC,merged$deseq2_log2FC,use='complete.obs')>0)

	return(merged)
}


count_de_genes=function(merged,M=c('up','down'),FDR=0.05){
	M=match.arg(M)
	if (M=='up'){
		merged%>%dplyr::filter(deseq2_log2FC>0,noiseq_log2FC>0,padj<FDR,prob>1-FDR)%>%nrow()
	} else {
		merged%>%dplyr::filter(deseq2_log2FC<0,noiseq_log2FC<0,padj<FDR,prob>1-FDR)%>%nrow()
	}
}


count_de_genes_at_many_FDRs=function(merged,FDR=c(0.001,0.01,0.05)){
	foreach(i=FDR,.combine='rbind')%do%{
		foreach(j=c('up','down'),.combine='rbind')%do%{
			n=count_de_genes(merged,j,i)
			data.table(n,direction=j,FDR=i)
		}
	}
}

plot_sig_genes=function(n_sig,add_label=FALSE,tissue_set=c("Artery - Coronary", 'Heart - Left Ventricle','Cells - Transformed fibroblasts')){
	n_sig[,FDR:=as.character(FDR)]
	p=ggplot(n_sig[tissue%in%tissue_set],
		aes(tissue,n,fill=direction,alpha=FDR,label=n))+
		geom_bar(stat='identity',position=position_dodge())+
		ylab('DE genes')+xlab('')+
		theme(axis.text.x=element_text(angle=45,hjust=1))+
		scale_fill_manual(values=c(down='blue',up='red'))+
		scale_alpha_manual(values=c(0.4,0.7,1.0))
	if (add_label) {
		p=p+geom_text(position=position_dodge(0.9),hjust=0,angle=90)
	}
	n_sig[,FDR:=as.numeric(FDR)]
	return(p)
}


# Main:
deseq2_list=read_deseq2_result(deseq2_fn)
noiseq_list=read_noiseq_result(noiseq_fn)

tissue = "Artery - Aorta"

n_sig=foreach(tissue=tissue_list,.combine='rbind')%dopar%{
	deseq2=select_tissue(deseq2_list,tissue)
	noiseq=select_tissue(noiseq_list,tissue)

	merged=merge_noiseq_and_deseq2(noiseq,deseq2,tissue)
	temp=count_de_genes_at_many_FDRs(merged,FDR=c(0.001,0.01,0.05))
	temp$tissue=tissue

	return(temp)
}

fwrite(n_sig,sprintf('%s/n_sig.txt',out_dir),sep='\t')

p1=plot_sig_genes(n_sig,add_label=FALSE)+
	scale_x_discrete(labels=c('Coronary artery','Fibroblast','Heart - left ventricle'))
p1b=plot_sig_genes(n_sig,add_label=FALSE,tissue_set=c("Artery - Coronary","Cells - Transformed fibroblasts"))+
	scale_x_discrete(labels=c('Coronary artery','Fibroblast'))
p2=plot_sig_genes(n_sig,add_label=FALSE,tissue_set=tissue_list)

pdf(sprintf('%s/n_sig.pdf',fig_dir))
p1;p1b;p2
dev.off()

saveRDS(list(p1=p1,p1b=p1b,p2=p2),sprintf('%s/fig.rda',out_dir))
