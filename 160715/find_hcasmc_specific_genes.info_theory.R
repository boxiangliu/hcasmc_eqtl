# Library:
library(cowplot)
library(data.table)
library(dplyr)

# Functions: 
subset_genes=function(x,keep){
	stopifnot(class(keep)=='character')
	y=x[Name%in%keep]
	return(y)
}

median2=function(x,coldata){
	w=data.frame()
	tissue=unique(coldata$tissue)
	for (i in tissue){
		soi=coldata[tissue==i,sample]
		y=x[,colnames(x)%in%soi,with=F]
		z=data.frame(apply(y,1,function(k) {median(k,na.rm=T)}))
		colnames(z)=i
		if (nrow(w)==0){
			w=z
		} else {
			w=cbind(w,z)
		}
	}
	return(w)
}

remove_zero_rows=function(x){
	y=x[rowSums(x)!=0,]
	return(y)
}

normalize=function(x){
	y=apply(x,1,sum)
	z=x/y
	return(z)
}

entropy=function(x){
	w=apply(x,1,function(y) {z=y[y!=0];-sum(z*log2(z))})
	return(w)
}

calculate_esi=function(x,tissue){
	# toi=tissue of interest
	# ot=other tissues
	stopifnot(class(x)=='data.frame')
	y=remove_zero_rows(x)
	toi=y[,tissue]
	ot=y[,colnames(y)!=tissue]
	toin=toi/rowSums(ot)
	otn=normalize(ot)
	e=entropy(otn)
	z=e-log2(toin)
	z=z[!is.infinite(z)]
	z=z[!is.na(z)]
	esi=1-z/max(z)
	return(esi)
}

calculate_esi_tissue_list=function(x,list_){
	esi_ls=list()
	stopifnot(all(list_%in%colnames(x)))
	for (i in list_){
		esi_ls[[i]]=calculate_esi(x,i)
	}
	return(esi_ls)
}

gid2gname=function(x,rowdata){
	x_id=names(x)
	x_name=rowdata$Description[match(x_id,rowdata$Name)]
	names(x)=x_name
	return(x)
}

list2df=function(x){
	y=list()
	for (i in names(x)){
		y[[i]]=data.frame(gene=names(x[[i]]),esi=x[[i]],tissue=i)
	}
	z=Reduce(rbind,y)
	w=dcast(z,gene~tissue,value.var='esi')
	return(w)
}

make_plot=function(rpkm,rowdata,coldata,gene_name,tissue_color){
	stopifnot(nrow(rpkm)==nrow(rowdata))
	stopifnot(ncol(rpkm)==nrow(coldata))
	x=data.frame(rpkm=t(rpkm[which(rowdata$Description==gene_name),]))
	x$sample=rownames(x)
	y=merge(x,coldata,by='sample')
	z=merge(y,tissue_color,by='tissue')
	z$log10_rpkm=log10(z$rpkm)
	z=z[!is.na(z$log10_rpkm),]
	z=z[is.finite(z$log10_rpkm),]
	filter=z%>%group_by(tissue)%>%summarise(num=n())
	z=z[z$tissue%in%filter[filter$num>=3,]$tissue,]
	w=unique(z[,c('tissue','color')])
	p=ggplot(z,aes(x=tissue,y=log10_rpkm,color=tissue,fill=tissue))+geom_violin()+scale_fill_manual(guide='none',breaks=w$tissue,values=w$color)+scale_color_manual(guide='none',breaks=w$tissue,values=w$color)+ylab('log10(RPKM)')+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+ggtitle(gene_name)
	return(p)
}

make_qqplot=function(esi_ls){
	tmpls=list()
	for (i in names(esi_ls)){
		z=qqplot(esi_ls[[i]],esi_ls[['HCASMC']],plot.it=F)
		tmpls[[i]]=data.frame(x=z$x,y=z$y,sample=i)
	}
	qqdata=Reduce(rbind,tmpls)
	p=ggplot(qqdata,aes(x=y,y=x,color=sample))+geom_line(alpha=ifelse(qqdata$sample=='HCASMC',1,0.2))+scale_color_manual(guide='none',breaks=tissue_color$tissue,values=tissue_color$color)+xlab('HCASMC quantiles')+ylab('GTEx quantiles')
	return(p)
}


# Main: 
# Read RPKM: 
x=fread('../processed_data/160715/combined.gtex.hcasmc.rpkm',header=T)
setnames(x,'9070202_Nextseq','9070202')


# Subset to protein coding gene and lncRNA:
keep=fread('../data/gtex/gencode.v19.genes.v6p.patched_contigs_genetypes.bed',header=F)%>%filter(V1%in%as.character(1:22))%>%filter(V6%in%c('protein_coding','lincRNA'))
x=subset_genes(x,keep$V5)


# Split x into rowdata, coldata and rpkm: 
rowdata=x[,.(Name,Description)]
rpkm=x[,c('Name','Description'):=NULL]
coldata=data.table(sample=colnames(rpkm))


# Read sample name to tissue table: 
sample_list=fread('../data/gtex/SAMPID_SMTSD.sorted.filtered.with_hcasmc.txt',header=F)
setnames(sample_list,c('sample','tissue'))
coldata=merge(coldata,sample_list)
coldata=coldata[match(colnames(rpkm),coldata$sample)]


# Calculate median of each tissue and each gene:
tissue_median=median2(rpkm,coldata)
rownames(tissue_median)=rowdata$Name


# Calculate ESI:
esi=calculate_esi(tissue_median,'HCASMC')
esi_gname=gid2gname(esi,rowdata)


# Save ESI: 
fwrite(data.table(gene_id=names(esi),esi=esi),'../processed_data/160715/esi.hcasmc.txt',sep='\t')


# Read tissue color: 
tissue_color=fread('shared/tissue_color.txt')[,c(1,5)]
colnames(tissue_color)=c('tissue','color')


# Plot some HCASMC-specific genes:
for (gene in c('MMP1','CXCL6','CALD1','VIM','MYH10','TPM4','RBP1')){
	p=make_plot(rpkm,rowdata,coldata,gene,tissue_color)
	save_plot(sprintf('../figures/160715/hcasmc_specific_gene.%s.pdf',gene),p)
}


# Calculate ESI for all tissues: 
esi_ls=calculate_esi_tissue_list(tissue_median,colnames(tissue_median))


# Save ESI for all tissues: 
esi_df=list2df(esi_ls)
fwrite(esi_df,'../processed_data/160715/esi.all_tissues.txt',sep='\t')


# Make qqplot:
p=make_qqplot(esi_ls)
save_plot('../figures/160715/esi_distribution_qqplot.pdf',p)


# Make RPKM vs rank plot: 
median_long=as.data.table(melt(tissue_median,variable.name='tissue',value.name='rpkm'))
median_long[,rank:=rank(-rpkm),by='tissue']
p=ggplot(median_long,aes(x=rank,y=log10(rpkm),color=tissue))+geom_line(alpha=ifelse(median_long$tissue=='HCASMC',1,0.2))+scale_color_manual(guide='none',breaks=tissue_color$tissue,values=tissue_color$color)
save_plot('../figures/160715/rpkm_vs_rank.pdf',p)