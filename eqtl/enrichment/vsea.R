# Perform variant set enrichment analysis
# Boxiang Liu (bliu2@stanford.edu)
# 2017-11-15

library(stringr)
library(data.table)
library(dplyr)
library(dtplyr)
library(cowplot)

source('/srv/persistent/bliu2/HCASMC_eQTL/scripts/VSEA/R/select_background_variants.R')
source('/srv/persistent/bliu2/HCASMC_eQTL/scripts/VSEA/R/select_LD_variants.R')
source('/srv/persistent/bliu2/HCASMC_eQTL/scripts/VSEA/R/variant_set_enrichment.R')

# Variables:
treeQTL_fn='../processed_data/rasqual/output_merged/treeQTL/eGenes.tsv'
in_dir='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/rasqual/output/'
vcf_dir='../processed_data/gwas_atacseq_overlap/prepare_vcf/'
anno_fn='../processed_data/sqtl/enrichment/snpEff/annotation/all.bed.gz'
out_dir='../processed_data/eqtl/enrichment/'
fig_dir='../figures/eqtl/enrichment/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

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

calc_odds_ratio=function(overlap_){
	message('INFO - calculating odds ratio')
	foreground_odds=overlap_[background_variant==foreground_variant,sum(loci_overlap)/.N]
	background_overlap=overlap_[background_variant!=foreground_variant]

	odds_ratio=foreach(i=1:500,.combine=c)%dopar%{
		bootstrap=background_overlap[sample(1:nrow(background_overlap),nrow(background_overlap),replace=TRUE)]
		background_odds=bootstrap[,sum(loci_overlap)/.N]
		foreground_odds/background_odds
	}

	mean=mean(odds_ratio)
	sd=sd(odds_ratio)
	return(data.table(mean=mean,sd=sd))
}

main=function(){
	# Select top variant:
	treeQTL=read_treeQTL(treeQTL_fn)
	top_eqtl=select_top_eqtl(treeQTL,1000,in_dir)

	snpid=unique(extract_snpid(top_eqtl))
	length(snpid) # 994

	# Select background variant:
	snpsnap=read_snpsnap()
	background_variants=select_background_variants(snpid,snpsnap,200)
	index_set=create_index_set(background_variants,snpsnap)
	fwrite(index_set,sprintf('%s/index_set.tsv',out_dir),sep='\t')

	# Select LD variant: 
	ld_set=select_LD_variants(index_set,out_dir,vcf_dir)
	fwrite(ld_set,sprintf('%s/ld_set.tsv',out_dir),sep='\t')

	# Calculate enrichment p-value:
	anno=read_annotation(anno_fn)
	enrichment=foreach(i=unique(anno$anno),.combine='rbind')%dopar%{
		message(i)
		overlap_=overlap(ld_set,anno[anno==i])
		pval=calc_enrichment(overlap_)
		odds_ratio=calc_odds_ratio(overlap_)
		data.table(anno=i,pval=pval,odds_ratio)
	}

	# Make plot:
	setorder(enrichment,-pval)
	enrichment[,anno:=factor(as.character(anno),level=as.character(anno))]
	p=ggplot(enrichment,aes(anno,-log10(pval)))+geom_point()+xlab('')+ylab('-log10(P-value)')+coord_flip()
	pdf(sprintf('%s/vsea.logp.pdf',fig_dir),height=4,width=4*1.3)
	p
	dev.off()

	setorder(enrichment,mean)
	enrichment[,anno:=factor(as.character(anno),level=as.character(anno))]
	p1=ggplot(enrichment,aes(anno,mean))+geom_pointrange(aes(ymin=mean-sd,ymax=mean+sd))+xlab('')+ylab('Odds ratio')+coord_flip()
	p2=ggplot(enrichment[!anno%in%c('splice_acceptor_variant','splice_donor_variant')],aes(anno,mean))+geom_pointrange(aes(ymin=mean-sd,ymax=mean+sd))+xlab('')+ylab('Odds ratio')+coord_flip()
	pdf(sprintf('%s/vsea.odds_ratio.pdf',fig_dir),height=4,width=4*1.3)
	p1;p2
	dev.off()
}

main()