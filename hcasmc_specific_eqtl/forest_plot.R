#!/usr/bin/r -i
library(cowplot)
library(meta)

# Variables:
tissue_color_file='shared/tissue_color.txt'

# Command args:
in_file=commandArgs(T)[1]
out_fig=commandArgs(T)[2]

# Read eqtl:
eqtl=read.table(in_file,header=F)
colnames(eqtl)=c('sid','pval','beta','se','tissue')

# Add tissue color:
tissue_color=read.table(tissue_color_file,header=T,sep='\t',comment.char='@')[,c('tissue_site_detail','tissue_site_detail_id','tissue_color_hex')]
colnames(tissue_color)=c('tissue_site_detail','tissue','color')
eqtl=merge(eqtl,tissue_color[,c('tissue','tissue_site_detail','color')],by='tissue')

# Calculate random effect size and se: 
re_result=metagen(eqtl$beta,eqtl$se,comb.random=T,comb.fixed=F)
re=with(re_result,data.frame(sid=eqtl$sid[1],pval=pval.random,beta=TE.random,se=seTE.random,tissue='RE Summary',tissue_site_detail='RE Summary',color='#000000'))
eqtl=rbind(eqtl,re)

# Make forest plot: 
eqtl$tissue_site_detail=factor(eqtl$tissue_site_detail,levels=rev(levels(eqtl$tissue_site_detail)))
p=ggplot(eqtl,aes(tissue_site_detail,beta))+
geom_errorbar(aes(ymax=beta+1.96*se,ymin=beta-1.96*se),color=eqtl$color,width=0.1,size=1)+
geom_point(size=1)+
geom_hline(yintercept=0,size=0.3)+
coord_flip()+
xlab('Tissue')+
ylab('Effect size')+
theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())

# Save plot:
save_plot(out_fig,p,base_height=6,base_width=4)

