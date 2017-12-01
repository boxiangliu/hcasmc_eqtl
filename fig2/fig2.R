# Make figure 2
# Boxiang Liu
# 2017-12-01

library(data.table)

# Variables:
in_dir='../processed_data/fig2/'
fig_dir='../figures/fig2/'
if (!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}

# Functions:
plot_enrichment=function(partition_heritability,title='',gtex_anno_fn='../data/gtex/gtex_tissue_colors.with_hcasmc.txt'){
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
		coord_flip()+
		ggtitle(title)
}

partition_heritability=fread(sprintf('%s/partition_heritability.txt',in_dir))
fig2a=plot_enrichment(partition_heritability)
save_plot(sprintf('%s/fig2.pdf',fig_dir),fig2a)