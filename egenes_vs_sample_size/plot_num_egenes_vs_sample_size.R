#!/usr/bin/env Rscript
# bosh liu
# durga
# compare the number of egenes discovered in hcasmc from those in gtex

# library
library(cowplot)
library(MASS)


# command args: 
gtex_file='../processed_data/160530/gtex.v6p.egenes.summary.txt'
hcasmc_file='../processed_data/160805/hcasmc.eqtl.pc4.peer8.perm.padj.txt'
figure_prefix='../figures/egenes_vs_sample_size/'


# read inputs:
gtex=read.table(gtex_file,sep='\t',header=T)
hcasmc=read.table(hcasmc_file,header=T)

# append hcasmc statistics onto gtex: 
eGenes=sum(hcasmc$qval<0.05)
Expressed=nrow(hcasmc)
Size=52
gtex=rbind(gtex,data.table(Tissue='HCASMC',eGenes,Size,Expressed,Pretty='HCASMC'))


# make scatter plot:
p=ggplot(gtex,aes(x=Size,y=eGenes,label=Tissue))+geom_point(size=2)+stat_smooth(formula=y~x,method='rlm')+xlab('Sample Size')+ylab('Number of eGenes (FDR < 0.05)')+geom_text(aes(label=ifelse(Tissue=='HCASMC','HCASMC','')),angle=90,hjust=0,vjust=0.5,nudge_y=200)
save_plot(paste0(figure_prefix,'num_egenes_vs_sample_size.fullsize.pdf'),p)


# make scatter plot for tissues with size smaller than 150:
gtex_small=gtex[Size<150,]
p2=ggplot(gtex_small,aes(x=Size,y=eGenes,label=Tissue))+geom_point(size=2)+stat_smooth(formula=y~x,method='rlm')+xlab('Sample Size')+ylab('Number of eGenes (FDR < 0.05)')+geom_text(aes(label=ifelse(Tissue=='HCASMC','HCASMC','')),angle=90,hjust=0,vjust=0.5,nudge_y=200)
save_plot(paste0(figure_prefix,'num_egenes_vs_sample_size.size150.pdf'),p2)
