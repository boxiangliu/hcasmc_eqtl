library(data.table)
library(stringr)
library(cowplot)
library(foreach)
library(doMC)
registerDoMC(15)

# Variables:
partition_heritability_dir='../processed_data/gwas_atacseq_overlap/ldscore_regression/partition_heritability_merged/'
fig_dir='../figures/gwas_atacseq_overlap/ldscore_regression/plot_heritability_merged/'
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
		geom_pointrange()+coord_flip()+ggtitle(title)
	return(p1)
}

plot_enrichment_pvalue=function(partition_heritability,title){
	setorder(partition_heritability,-Enrichment_p)
	partition_heritability[,tissue:=factor(tissue,tissue)]
	p1=ggplot(partition_heritability,aes(x=tissue,y=-log10(Enrichment_p)))+
		geom_bar(stat='identity')+coord_flip()+ggtitle(title)+ylab('-log10(Enrichment P-value)')
	return(p1)
}

tissue_set_list=c('all_tissue','jaccard_similarity_0.3','jaccard_similarity_0.4','jaccard_similarity_0.5')
pdf(sprintf('%s/heritability_nobaseline.pdf',fig_dir),height=40)

for (tissue_set in tissue_set_list){
	partition_heritability=read_heritability(partition_heritability_dir,tissue_set)
	p=plot_enrichment_pvalue(partition_heritability,sprintf('%s',tissue_set))
	print(p)
}
dev.off()



partition_heritability[,Enrichment2:=Prop._h2/Prop._SNPs]
ggplot(partition_heritability,aes(Enrichment,Enrichment2))+geom_point()

partition_heritability[,Enrichment_z:=abs(Enrichment)/Enrichment_std_error]
ggplot(partition_heritability,aes(Enrichment_z,-log10(Enrichment_p)))+geom_point()

ggplot(partition_heritability,aes(Enrichment,1/Prop._SNPs))+geom_point()
ggplot(partition_heritability,aes(Enrichment,Prop._h2))+geom_point()
ggplot(partition_heritability,aes(-log10(Enrichment_p),Prop._h2))+geom_point()
ggplot(partition_heritability,aes(-log10(Enrichment_p),Prop._SNPs))+geom_point()
ggplot(partition_heritability,aes(Enrichment_z,`Coefficient_z-score`))+geom_point()
ggplot(partition_heritability,aes(-log10(Enrichment_p),-log10(Coefficient_p)))+geom_point()