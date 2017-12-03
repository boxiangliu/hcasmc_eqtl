# Perform variant set enrichment analysis
# Boxiang Liu (bliu2@stanford.edu)
# 2017-11-15

library(stringr)
library(data.table)
library(dplyr)
library(dtplyr)
library(cowplot)
library(vsea)
library(ggrepel)
library(foreach)
library(doMC)
registerDoMC(10)

# Variables:
in_fn='../processed_data/gwas_atacseq_overlap/tmp/ld_set.tsv'
anno_dir='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_released_adult_tissue_group/'
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

make_enrichment_plot=function(enrichment){
	# Make plot:
	gtex_anno_fn='../data/gtex/gtex_tissue_colors.with_hcasmc.txt'
	gtex_anno=fread(gtex_anno_fn)
	metadata_fn='../data/encode/dnase_seq/metadata.tsv'
	metadata=fread(metadata_fn)

	x=metadata[`Biosample type`%in%c('tissue','primary cell')&`Biosample life stage`%in%c('adult'),list(`GTEx tissue`,`Tissue group`)]
	y=unique(merge(x,gtex_anno,by.x='GTEx tissue',by.y='tissue_site_detail',all.x=TRUE,all.y=FALSE)[,list(`Tissue group`,`tissue_color_hex`)])
	z=rbind(y,data.table(`Tissue group`='HCASMC',tissue_color_hex='#FF0066'))
	color=z$tissue_color_hex
	names(color)=z$`Tissue group`
	enrichment[,label:=ifelse(rank(pval)<=5|rank(-mean)<=5,anno,'')]

	ggplot(enrichment,aes(-log10(pval),mean,label=label,color=anno))+
		geom_point(size=2)+geom_text_repel(color='black',segment.alpha=0,nudge_x=ifelse(enrichment$label=='Artery Endothelial Cell',-0.7,0))+
		xlab('-log10(P-value)')+ylab('Odds ratio')+
		scale_color_manual(values=color,guide='none')
}

main=function(in_fn,anno_dir,fig_dir,fig_prefix,out_path){
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
		n_overlap=count_overlap(overlap_)
		data.table(anno=i,pval=pval,odds_ratio,n_overlap)
	}
	if (!is.null(out_path)) fwrite(enrichment,out_path,sep='\t')

	p0=make_enrichment_plot(enrichment)
	save_plot(sprintf('%s/%s.enrichment.pdf',fig_dir,fig_prefix),p0)

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

	return(enrichment)
}

main(in_fn,'../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_released_adult/',fig_dir,'vsea.merge_peaks_released_adult',NULL)
main(in_fn,'../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_released_adult_tissue_group/',fig_dir,'vsea.merge_peaks_released_adult_tissue_group',sprintf('%s/merge_peaks_released_adult_tissue_group.enrichment.txt',out_dir))