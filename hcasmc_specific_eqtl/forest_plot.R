#!/usr/bin/r -i
library(cowplot)
library(meta)


in_file=commandArgs(T)[1]
out_fig=commandArgs(T)[2]
eqtl=read.table(in_file,header=F)
colnames(eqtl)=c('sid','pval','beta','se','tissue')
re_result=metagen(eqtl$beta,eqtl$se,comb.random=T,comb.fixed=F)
re=with(re_result,data.frame(sid=eqtl$sid[1],pval=pval.random,beta=TE.random,se=seTE.random,tissue='RE Summary'))
eqtl=rbind(eqtl,re)
eqtl$tissue=factor(eqtl$tissue,levels=rev(levels(eqtl$tissue)))
p=ggplot(eqtl,aes(rev(tissue),beta))+
geom_point(size=1)+
geom_hline(yintercept=0,size=0.3)+
geom_errorbar(aes(ymax=beta+1.96*se,ymin=beta-1.96*se),width=0.1)+
coord_flip()+
xlab('Tissue')+
ylab('Effect size')

save_plot(out_fig,p)
