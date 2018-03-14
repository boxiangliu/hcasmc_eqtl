library(stringr)
library(data.table)
library(dplyr)
library(dtplyr)
library(cowplot)
library(foreach)
library(doMC)
registerDoMC(50)

treeQTL_fn='../processed_data/rasqual/output_merged/treeQTL/eGenes_level0.05.tsv'
anno_fn='../processed_data/sqtl/enrichment/snpEff/annotation/all.bed.gz'

read_treeQTL=function(in_fn){
	if (str_detect(in_fn,'.gz')){
		x=fread(sprintf('zcat %s',in_fn))
	} else {
		x=fread(in_fn)
	}
	return(x)
}

read_rasqual=function(fn){
	if (str_detect(fn,'.gz')){
		x=fread(sprintf('zcat %s',fn),select=c(1:7,11:12,25),col.names=c('fid','sid','chr','pos','ref','alt','af','chisq','pi','r2_rSNP'))
	} else {
		x=fread(fn,select=c(1:7,11:12,25),col.names=c('fid','sid','chr','pos','ref','alt','af','chisq','pi','r2_rSNP'))
	}
	return(x)
}

select_top_eqtl=function(treeQTL,n,in_dir){
	if (n>nrow(treeQTL)){
		warning('n is larger than nrow(treeQTL); using all eQTLs.')
		n=nrow(treeQTL)
	}
	x=treeQTL[rank(fam_p)<=n,]
	top_eqtl=foreach(i=x$family,.combine='rbind')%dopar%{
		message('INFO - ',i)
		fn=list.files('../processed_data/rasqual/output/',pattern=i,full.names=TRUE,recursive=TRUE)
		stopifnot(length(fn)==1)
		eqtl=read_rasqual(fn)
		eqtl=eqtl[r2_rSNP>0.8]
		eqtl[,max_chisq:=max(chisq)]
		y=eqtl[chisq==max_chisq,]
		if (nrow(y)>1){
			set.seed(42)
			y=y[sample(nrow(y),1)]
		}
		return(y)
	}
	return(top_eqtl)
}

extract_snpid=function(x){
	y=x[,list(chr,pos)]
	y[,chr:=str_replace(chr,'chr','')]
	return(y[,paste(chr,pos,sep=':')])
}

read_annotation=function(anno_fn){
	if (str_detect(anno_fn,'gz')){
		x=fread(sprintf('zcat %s',anno_fn))
	} else {
		x=fread(anno_fn)
	}
	setnames(x,c('chr','start','end','anno'))
	if (!str_detect(x$chr[1],'chr')){
		x[,chr:=paste0('chr',chr)]
	}
	x[,snpID:=paste(chr,start,sep=':')]
	return(x)
}

treeQTL=read_treeQTL(treeQTL_fn)
top_eqtl=select_top_eqtl(treeQTL,nrow(treeQTL),in_dir)
snpid=unique(extract_snpid(top_eqtl))
anno=read_annotation(anno_fn)
