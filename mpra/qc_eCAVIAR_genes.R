#!/usr/env Rscript 
library(data.table)
library(dplyr)
library(dtplyr)
library(cowplot)
args=commandArgs(T)
in_file=args[1]
fig_file=args[2]
# in_file='../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_genes.txt'
# fig_file='../figures/mpra/eCAVIAR/qc_eCAVIAR_colocalized_genes.pdf'

genes=fread(in_file)
clpp_cutoff=0.01
genes[,colocalized:=clppSum>clpp_cutoff]
num_tested=nrow(genes)
num_colocalized=sum(genes$clppSum>clpp_cutoff)

# Plot distribution of CLPP scores:
p1=ggplot(genes,aes(x=clppSum))+geom_histogram(binwidth=0.005)+xlab('Sum of CLPP score')+geom_vline(xintercept=0.01,color='red',linetype=2)+scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0,1.2))+annotate(geom='text',x=0.4,y=40,label=paste0(num_colocalized,'/',num_tested,' genes with sum of CLPP > 0.01'))


# Plot distribution of putatively colocalized variants: 
p2=ggplot(genes,aes(x=gwasSet,color=colocalized))+geom_density()+xlab('GWAS')
p3=ggplot(genes,aes(x=eqtlSet,color=colocalized))+geom_density()+xlab('eQTL')
p4=ggplot(genes,aes(x=intersectSet,fill=colocalized))+geom_histogram(position='dodge',binwidth=1)+xlab('Intersection')
p5=ggplot(genes,aes(x=unionSet,color=colocalized))+geom_density()+xlab('Union')
p6=plot_grid(p2,p3,p4,p5,labels=LETTERS[1:4])


# save plot:
pdf(fig_file)
print(p1)
print(p6)
dev.off()
