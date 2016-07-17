#!/usr/bin/env Rscript
# bosh liu
# durga
# compare the number of egenes discovered in hcasmc from those in gtex

# library
library(cowplot)


# command args: 
args=commandArgs(T)
gtex_file=args[1]
hcasmc_file=args[2]
figure=args[3]


# read inputs:
gtex=read.table(gtex_file,sep='\t',header=T)
hcasmc=read.table(hcasmc_file,header=T)

# append hcasmc statistics onto gtex: 
eGenes=sum(hcasmc$qvalues<0.05)
Expressed=nrow(hcasmc)
Size=52
gtex=rbind(gtex,data.table(Tissue='HCASMC',eGenes,Size,Expressed,Pretty='HCASMC'))


# make scatter plot:
p=ggplot(gtex,aes(x=Size,y=eGenes,label=Tissue))+geom_point(size=2)+geom_smooth(formula=y~x)+xlab('Sample Size')+ylab('Number of eGenes (FDR < 0.05)')+geom_text(aes(label=ifelse(Tissue=='HCASMC','HCASMC','')),angle=90,hjust=0,vjust=0.5,nudge_y=200)
save_plot(figure,p)