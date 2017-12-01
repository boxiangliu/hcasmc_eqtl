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
in_fn='../processed_data/gwas_atacseq_overlap/tmp/ld_set.tsv'
anno_dir='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_released_adult/'
out_dir='../processed_data/gwas_atacseq_overlap/vsea/'
fig_dir='../figures/gwas_atacseq_overlap/vsea/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

read_ld_set=function(in_fn){
	ld_set=fread(in_fn)
	ld_set[,c('background_variant','foreground_variant'):=list(loci_index,gwas_index)]
	ld_set[,c('start','end'):=list(pos,pos)]
	return(ld_set)
}

read_annotation=function(in_dir){
	fn=list.files(in_dir,pattern='bed')
	x=foreach(f=fn,.combine='rbind')%dopar%{
		sample=str_replace(f,'.merged.bed','')
		print(sprintf('INFO - %s',sample))
		dhs=fread(sprintf('%s/%s',in_dir,f))
		dhs$anno=sample
		return(dhs)
	}
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

count_overlap=function(in_dir,ld_set){
	fn=list.files(in_dir,pattern='bed')
	gwas_set=ld_set[snpID==gwas_index]
	setkey(gwas_set,chr,start,end)
	n_overlap=foreach(f=fn,.combine='rbind')%dopar%{
		sample=str_replace(f,'.merged.bed','')
		print(sprintf('INFO - %s',sample))
		dhs=fread(sprintf('%s/%s',in_dir,f))
		setkey(dhs,chr,start,end)


		# Overlap: 
		overlap=unique(foverlaps(gwas_set,dhs[,list(chr,start,end)]))
		overlap[,c('i.start','i.end'):=NULL]


		overlap[,snp_overlap:=!is.na(start)]
		overlap[,loci_overlap:=any(snp_overlap),by='loci_index']
		overlap[,c('start','end'):=NULL]
		overlap=unique(overlap)
		stopifnot(nrow(overlap)==nrow(gwas_set))


		overlap=overlap[ld_proxy==FALSE,]
		n_overlap=unlist(overlap[,list(n=sum(loci_overlap))])
		data.table(sample,n_overlap)
	}
	setorder(n_overlap,-n_overlap)
}


main=function(in_fn,anno_dir,fig_dir,fig_prefix){
	message('LD set: ', in_fn)
	message('annotation: ', anno_dir)

	# Read LD variant: 
	ld_set=read_ld_set(in_fn)

	# Calculate enrichment p-value:
	anno=read_annotation(anno_dir)
	enrichment=foreach(i=unique(anno$anno),.combine='rbind')%dopar%{
		message(i)
		overlap_=overlap(ld_set,anno[anno==i])
		pval=calc_enrichment(overlap_)
		odds_ratio=calc_odds_ratio(overlap_)
		data.table(anno=i,pval=pval,odds_ratio)
	}

	# n_overlap=count_overlap(anno_dir,ld_set)
	# enrichment[anno%in%n_overlap[n_overlap>=1,sample]]

	# Make plot:
	setorder(enrichment,-pval)

	enrichment[,anno:=factor(as.character(anno),level=as.character(anno))]
	p=ggplot(enrichment,aes(anno,-log10(pval)))+geom_point()+xlab('')+ylab('-log10(P-value)')+coord_flip()
	pdf(sprintf('%s/%s.logp.pdf',fig_dir,fig_prefix),height=8,width=8*1.3)
	print(p)
	dev.off()


	setorder(enrichment,mean)
	enrichment[,anno:=factor(as.character(anno),level=as.character(anno))]
	p1=ggplot(enrichment,aes(anno,mean))+geom_pointrange(aes(ymin=mean-sd,ymax=mean+sd))+xlab('')+ylab('Odds ratio')+coord_flip()
	pdf(sprintf('%s/%s.odds_ratio.pdf',fig_dir,fig_prefix),height=8,width=8*1.3)
	print(p1)
	dev.off()
}

main(in_fn,'../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_released_adult/',fig_dir,'vsea.merge_peaks_released_adult')
main(in_fn,'../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_released_adult_tissue_group/',fig_dir,'vsea.merge_peaks_released_adult_tissue_group')