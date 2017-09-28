library(data.table)
library(stringr)
library(cowplot)
library(foreach)
library(doMC)
registerDoMC(15)

# Variables:
partition_heritability_fn='../processed_data/gwas_gene_overlap/ldscore_regression/partition_heritability_merged/cad.merged.results'
partition_heritability_nobaseline_fn='../processed_data/gwas_gene_overlap/ldscore_regression/partition_heritability_merged/cad.merged.nobaseline.results'
fig_dir='../figures/gwas_gene_overlap/ldscore_regression/plot_heritability_merged/'
if (!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}

#------------- Partition heritability ----------#
# Read partitioned heritability: 
partition_heritability=fread(partition_heritability_fn)[1:22]
partition_heritability[,tissue:=str_replace_all(str_replace(Category,'L2_0',''),'_',' ')]
partition_heritability[,Coefficient_p:=2*pnorm(-abs(`Coefficient_z-score`))]


# Plot partitioned heritability for each tissue:
setorder(partition_heritability,Coefficient_p)
partition_heritability[,tissue:=factor(tissue,tissue)]
p1=ggplot(partition_heritability,aes(x=tissue,y=-log10(Coefficient_p)))+
	geom_point()+coord_flip()

setorder(partition_heritability,Enrichment)
partition_heritability[,tissue:=factor(tissue,tissue)]
p2=ggplot(partition_heritability,aes(x=tissue,y=Enrichment,ymin=Enrichment-Enrichment_std_error,ymax=Enrichment+Enrichment_std_error))+
	geom_pointrange()+coord_flip()

pdf(sprintf('%s/heritability.pdf',fig_dir))
p1;p2
dev.off()

#-------------- Partition heritability with baseline ------------#
# Read partitioned heritability: 
partition_heritability=fread(partition_heritability_nobaseline_fn)[1:22]
partition_heritability[,tissue:=str_replace_all(str_replace(Category,'L2_0',''),'_',' ')]
partition_heritability[,Coefficient_p:=2*pnorm(-abs(`Coefficient_z-score`))]


# Plot partitioned heritability for each tissue:
setorder(partition_heritability,Coefficient_p)
partition_heritability[,tissue:=factor(tissue,tissue)]
p3=ggplot(partition_heritability,aes(x=tissue,y=-log10(Coefficient_p)))+
	geom_point()+coord_flip()

setorder(partition_heritability,Enrichment)
partition_heritability[,tissue:=factor(tissue,tissue)]
p4=ggplot(partition_heritability,aes(x=tissue,y=Enrichment,ymin=Enrichment-Enrichment_std_error,ymax=Enrichment+Enrichment_std_error))+
	geom_pointrange()+coord_flip()

pdf(sprintf('%s/heritability_nobaseline.pdf',fig_dir))
p3;p4
dev.off()