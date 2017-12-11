# Make locuscatter plots
# Boxiang Liu
# 2017-12-07
library(R.utils)
library(data.table)
library(cowplot)
library(ggrepel)

args=commandArgs(T,T,
	defaults=list(marker_col1='rsid',
				  pval_col1='pval',
				  marker_col2='rsid',
				  pval_col2='pval',
				  title1='Study 1',
				  title2='Study 2',
				  snp=NULL))

in_fn1=args$f1
marker_col1=args$marker_col1
pval_col1=args$pval_col1
title1=args$title1

in_fn2=args$f2
marker_col2=args$marker_col2
pval_col2=args$pval_col2
title2=args$title2
snp=args$snp

read_metal=function(in_fn,marker_col,pval_col){
	if (grepl('.gz',in_fn)){
		d=fread(sprintf('gunzip -c %s',in_fn))
	} else {
		d=fread(in_fn)
	}
	setnames(d,c(marker_col,pval_col),c('rsid','pval'))
	d[,list(rsid,pval,logp=-log10(pval))]
}


make_locuscatter=function(merged,title1,title2,snp){
	if (is.null(snp)){
		snp=merged[pval1==min(pval1),rsid]
	} else {
		if(!snp%in%merged$rsid){
			stop(sprintf('%s not found in %s',snp,in_fn1))
		}
	}
	color=ifelse(merged$rsid==snp,'purple','black')
	names(color)=merged$rsid
	shape=ifelse(merged$rsid==snp,18,16)
	names(shape)=merged$rsid
	size=ifelse(merged$rsid==snp,5,2)
	names(size)=merged$rsid
	merged[,label:=ifelse(rsid==snp,rsid,'')]

	p=ggplot(merged,aes(logp1,logp2))+
		geom_point(aes(color=rsid,size=rsid,shape=rsid),alpha=0.8)+
		xlab(title1)+ylab(title2)+
		scale_color_manual(values=color,guide='none')+
		scale_shape_manual(values=shape,guide='none')+
		scale_size_manual(values=size,guide='none')+
		geom_text_repel(aes(label=label))
	return(p)
}

in_fn1='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/rasqual/output_pval/chr15/ENSG00000182511.7_FES.pval.txt'
marker_col1='rsid'
pval_col1='pval'

in_fn2='/srv/persistent/bliu2/HCASMC_eQTL/data//gwas/ukbb/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz'
marker_col2='snptestid'
pval_col2='p-value_gc'

d1=read_metal(in_fn1,marker_col1,pval_col1)
d2=read_metal(in_fn2,marker_col2,pval_col2)
merged=merge(d1,d2,by='rsid',suffixes=c('1','2'),all=FALSE)

p=make_locuscatter(merged,'eQTL','UKBB',NULL)
pdf('1.pdf',height=4,width=4)
print(p)
dev.off()