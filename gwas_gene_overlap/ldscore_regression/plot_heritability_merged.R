library(data.table)
library(stringr)
library(cowplot)
library(foreach)
library(doMC)
registerDoMC(15)

# Variables:
partition_heritability_dir='../processed_data/gwas_gene_overlap/ldscore_regression/partition_heritability_merged/'
fig_dir='../figures/gwas_gene_overlap/ldscore_regression/plot_heritability_merged/'
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
plot_enrichment=function(partition_heritability,title){
	setorder(partition_heritability,Enrichment)
	partition_heritability[,tissue:=factor(tissue,tissue)]
	p1=ggplot(partition_heritability,aes(x=tissue,y=Enrichment,ymin=Enrichment-Enrichment_std_error,ymax=Enrichment+Enrichment_std_error))+
		geom_pointrange()+xlab('')+coord_flip()+ggtitle(title)
	return(p1)
}

# Plot partitioned heritability for each tissue for the supplement:
plot_enrichment_for_supplement=function(partition_heritability,title){
	setorder(partition_heritability,Enrichment)
	partition_heritability[,tissue:=factor(tissue,unique(tissue))]
	ggplot(partition_heritability,aes(x=tissue,y=Enrichment,ymin=Enrichment-Enrichment_std_error,ymax=Enrichment+Enrichment_std_error,color=gene_set))+
		geom_pointrange(position=position_dodge(width=0.7))+xlab('')+coord_flip()+ggtitle(title)
}

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


# Plot for supplementary material:
partition_heritability=foreach(i=c('top1000','top2000','top4000'),.combine='rbind')%do%{
	temp=read_heritability(partition_heritability_dir,sprintf('%s/%s','tissue_specific_gene_no_sm',i))
	data.table(temp,gene_set=i)
}

p2=plot_enrichment_for_supplement(partition_heritability,'')
pdf(sprintf('%s/heritability_nobaseline.supplement.pdf',fig_dir))
p2
dev.off()

