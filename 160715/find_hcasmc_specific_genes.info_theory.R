# Library:
library(cowplot)
library(data.table)

# Functions: 
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
	message(min(e,na.rm=T))
	message(min(-log2(toin),na.rm=T))
	z=e-log2(toin)
	z=z[!is.infinite(z)]
	z=z[!is.na(z)]
	message(max(z))
	message(min(z))
	esi=1-z/max(z)
	return(esi)
}


make_plot=function(rpkm,rowdata,coldata,gene_id,tissue_color){
	x=data.frame(rpkm=t(rpkm[which(rowdata$Name==gene_id),]))
	x$sample=rownames(x)
	y=merge(x,coldata,by='sample')
	z=merge(y,tissue_color,by='tissue')
	w=unique(z[,c('tissue','color')])
	p=ggplot(z,aes(x=tissue,y=log10(rpkm),color=tissue,fill=tissue))+geom_violin()+geom_boxplot(width=0.1,outlier.shape=NA)+scale_fill_manual(guide='none',breaks=w$tissue,values=w$color)+scale_color_manual(guide='none',breaks=w$tissue,values=w$color)+ylab('log10(RPKM)')+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
	return(p)
}

# Main: 
# Read RPKM: 
x=fread('../processed_data/160715/combined.gtex.hcasmc.rpkm',header=T)
setnames(x,'9070202_Nextseq','9070202')
rowdata=rpkm[,.(Name,Description)]
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


# Read tissue color: 
tissue_color=fread('shared/tissue_color.txt')[,c(1,5)]
colnames(tissue_color)=c('tissue','color')


# Plot:
p1=make_plot(rpkm,rowdata,coldata,'ENSG00000238172.1',tissue_color)
p2=make_plot(rpkm,rowdata,coldata,'ENSG00000213772.3',tissue_color)
save_plot('../figures/160715/hcasmc_specific_gene.ENSG00000238172.1.pdf',p1)
save_plot('../figures/160715/hcasmc_specific_gene.ENSG00000213772.3.pdf',p2)