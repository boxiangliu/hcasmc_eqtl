# Make figure 2
# Boxiang Liu
# 2017-12-01

library(data.table)
library(cowplot)
library(ggrepel)

# Variables:
in_dir='../processed_data/fig2/'
fig_dir='../figures/fig2/'
if (!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}

# Functions:
plot_heritability=function(partition_heritability,title='',gtex_anno_fn='../data/gtex/gtex_tissue_colors.with_hcasmc.txt'){
	gtex_anno=fread(gtex_anno_fn)
	color=gtex_anno$tissue_color_hex
	names(color)=gtex_anno$abbreviation
	partition_heritability=merge(partition_heritability,gtex_anno[,list(tissue_site_detail,tissue_color_hex,abbreviation)],by.x='tissue',by.y='tissue_site_detail',all.x=TRUE,all.y=FALSE)
	setorder(partition_heritability,Enrichment)
	partition_heritability[,abbreviation:=factor(abbreviation,unique(abbreviation))]
	ggplot(partition_heritability,aes(x=abbreviation,y=Enrichment,
		ymin=Enrichment-Enrichment_std_error,ymax=Enrichment+Enrichment_std_error,color=abbreviation))+
		geom_linerange(color='grey',size=1)+
		geom_point(size=2)+xlab('')+
		scale_color_manual(values=color,guide='none')+
		theme(axis.text.y=element_text(color=ifelse(partition_heritability$abbreviation=='HCASMC','purple','black')))+
		coord_flip()
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
		geom_point(size=2)+geom_text_repel(color=ifelse(enrichment$label=='HCASMC','purple','black'),segment.alpha=0,nudge_x=ifelse(enrichment$label=='Artery Endothelial Cell',-0.7,0))+
		xlab('-log10(P-value)')+ylab('Odds ratio')+
		scale_color_manual(values=color,guide='none')
}


partition_heritability=fread(sprintf('%s/partition_heritability.txt',in_dir))
enrichment=fread('../processed_data/gwas_atacseq_overlap/vsea/merge_peaks_released_adult_tissue_group.enrichment.txt')

fig2a=plot_heritability(partition_heritability)
fig2b=make_enrichment_plot(enrichment)
fig2=plot_grid(fig2a,fig2b, labels = c("A", "B"))

save_plot(sprintf('%s/fig2.pdf',fig_dir),fig2,base_width=8)