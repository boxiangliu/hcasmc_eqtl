library(data.table)
library(stringr)
library(cowplot)
library(foreach)
library(doMC)
registerDoMC(15)

# Variables:
partition_heritability_dir='../processed_data/gwas_gene_overlap/ldscore_regression/partition_heritability_merged_scz/'
fig_dir='../figures/gwas_gene_overlap/ldscore_regression/partition_heritability_merged_scz/'
tissue_set = 'tissue_specific_gene'
experiment = 'top2000'

if (!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}

# Functions: 

# Read partitioned heritability: 
read_heritability=function(partition_heritability_dir,experiment){
	partition_heritability_fn=sprintf('%s/%s/cad.merged.nobaseline.results',partition_heritability_dir,experiment)
	partition_heritability=fread(partition_heritability_fn)
	partition_heritability[,tissue:=str_replace_all(str_replace(Category,'L2_0',''),'_',' ')]
	partition_heritability[,Coefficient_p:=2*pnorm(-abs(`Coefficient_z-score`))]
	return(partition_heritability)
}

# Plot partitioned heritability for each tissue:
plot_enrichment=function(partition_heritability,title='',gtex_anno_fn='../data/gtex/gtex_tissue_colors.with_hcasmc.txt'){
	gtex_anno=fread(gtex_anno_fn)
	color=gtex_anno$tissue_color_hex
	names(color)=gtex_anno$abbreviation
	partition_heritability=merge(partition_heritability,gtex_anno[,list(tissue_site_detail,tissue_color_hex,abbreviation)],by.x='tissue',by.y='tissue_site_detail',all.x=TRUE,all.y=FALSE)
	setorder(partition_heritability,Enrichment)
	partition_heritability[,abbreviation:=factor(abbreviation,unique(abbreviation))]
	ggplot(partition_heritability,aes(x=abbreviation,y=Enrichment,
		ymin=Enrichment-Enrichment_std_error,ymax=Enrichment+Enrichment_std_error,color=abbreviation))+
		geom_linerange(color='grey',size=1.5)+
		geom_point(size=3)+xlab('')+
		scale_color_manual(values=color,guide='none')+
		theme(axis.text.y=element_text(color=ifelse(partition_heritability$abbreviation=='HCASMC','purple','black')))+
		coord_flip()+
		ggtitle(title)
}

# Plot coefficient: 
plot_coefficient=function(partition_heritability,title='',gtex_anno_fn='../data/gtex/gtex_tissue_colors.with_hcasmc.txt'){
	gtex_anno=fread(gtex_anno_fn)
	color=gtex_anno$tissue_color_hex
	names(color)=gtex_anno$abbreviation
	partition_heritability=merge(partition_heritability,gtex_anno[,list(tissue_site_detail,tissue_color_hex,abbreviation)],by.x='tissue',by.y='tissue_site_detail',all.x=TRUE,all.y=FALSE)
	setorder(partition_heritability,`Coefficient_z-score`)
	partition_heritability[,abbreviation:=factor(abbreviation,unique(abbreviation))]
	ggplot(partition_heritability,aes(x=abbreviation,y=`Coefficient_z-score`,color=abbreviation))+
		geom_point(size=3)+xlab('')+
		scale_color_manual(values=color,guide='none')+
		theme(axis.text.y=element_text(color=ifelse(partition_heritability$abbreviation=='HCASMC','purple','black')))+
		coord_flip()+
		ggtitle(title)
}

partition_heritability=read_heritability(partition_heritability_dir,sprintf('%s/%s',tissue_set,experiment))
p = plot_enrichment(partition_heritability,sprintf('%s\n%s',tissue_set,experiment))
p2 = plot_coefficient(partition_heritability,sprintf('%s\n%s',tissue_set,experiment))
pdf(sprintf('%s/heritability_nobaseline.pdf',fig_dir))
print(p);print(p2)
dev.off()