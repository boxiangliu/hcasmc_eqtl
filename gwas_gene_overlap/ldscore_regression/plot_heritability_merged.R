library(data.table)
library(stringr)
library(cowplot)
library(foreach)
library(doMC)
registerDoMC(15)

# Variables:
partition_heritability_dir='../processed_data/gwas_gene_overlap/ldscore_regression/partition_heritability_merged/'
fig_dir='../figures/gwas_gene_overlap/ldscore_regression/plot_heritability_merged/'
ms_dir='../processed_data/fig2/'
if (!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}
if (!dir.exists(ms_dir)){dir.create(ms_dir,recursive=TRUE)}

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

# Plot partitioned heritability for each tissue for the supplement:
plot_enrichment_for_supplement=function(partition_heritability,title,gtex_anno_fn='../data/gtex/gtex_tissue_colors.with_hcasmc.txt'){
	gtex_anno=fread(gtex_anno_fn)
	partition_heritability=merge(partition_heritability,gtex_anno[,list(tissue_site_detail,tissue_color_hex,abbreviation)],by.x='tissue',by.y='tissue_site_detail',all.x=TRUE,all.y=FALSE)

	setorder(partition_heritability,Enrichment)
	partition_heritability[,tissue:=factor(tissue,unique(tissue))]
	partition_heritability[,abbreviation:=factor(abbreviation,unique(abbreviation))]
	ggplot(partition_heritability,aes(x=abbreviation,y=Enrichment,
		ymin=Enrichment-Enrichment_std_error,ymax=Enrichment+Enrichment_std_error,color=gene_set,group=gene_set))+
		geom_linerange(size=1.5,position=position_dodge(width=0.7))+
		geom_point(size=3,position=position_dodge(width=0.7))+
		xlab('')+coord_flip()+ggtitle(title)+
		scale_color_discrete(name='Gene Set',labels=c('Top 1000','Top 2000','Top 4000'))
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

#------------#
# enrichment #
#------------#
# Plot all combinations:
tissue_set_list=c('tissue_specific_gene','tissue_specific_gene_no_sm','tissue_specific_gene_no_sm_no_blood')
experiment_list=c('top200','top500','top1000','top2000','top4000','top8000')
pdf(sprintf('%s/heritability_nobaseline.pdf',fig_dir))
for (tissue_set in tissue_set_list){
	p=foreach(experiment=experiment_list,.final=function(x) setNames(x,experiment_list))%dopar%{
		partition_heritability=read_heritability(partition_heritability_dir,sprintf('%s/%s',tissue_set,experiment))
		plot_enrichment(partition_heritability,sprintf('%s\n%s',tissue_set,experiment))
	}
	print(p)
}
dev.off()


# Save minimum reproducible example for Fig. 2:
partition_heritability=read_heritability(partition_heritability_dir,'tissue_specific_gene_no_sm/top2000')
fwrite(partition_heritability,sprintf('%s/partition_heritability.txt',ms_dir),sep='\t')


# Plot for supplementary material:
partition_heritability=foreach(i=c('top1000','top2000','top4000'),.combine='rbind')%do%{
	temp=read_heritability(partition_heritability_dir,sprintf('%s/%s','tissue_specific_gene_no_sm',i))
	data.table(temp,gene_set=i)
}

p2=plot_enrichment_for_supplement(partition_heritability,'')
pdf(sprintf('%s/heritability_nobaseline.supplement.pdf',fig_dir))
p2
dev.off()

#---------#
# z-score #
#---------#

# Plot all combinations:
tissue_set_list=c('tissue_specific_gene','tissue_specific_gene_no_sm','tissue_specific_gene_no_sm_no_blood')
experiment_list=c('top200','top500','top1000','top2000','top4000','top8000')
pdf(sprintf('%s/heritability_nobaseline_zscore.pdf',fig_dir))
for (tissue_set in tissue_set_list){
	p=foreach(experiment=experiment_list,.final=function(x) setNames(x,experiment_list))%dopar%{
		partition_heritability=read_heritability(partition_heritability_dir,sprintf('%s/%s',tissue_set,experiment))
		plot_coefficient(partition_heritability,sprintf('%s\n%s',tissue_set,experiment))
	}
	print(p)
}
dev.off()