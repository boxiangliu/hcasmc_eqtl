#!/usr/bin/env Rscript
# bosh liu
# durga
# plot the number of egenes vs number of genotype PCs and number of PEER factors


# library
library(stringr)
library(gtools)
library(cowplot)


# command args: 
args=commandArgs(T)
input_file=args[1]
figure_path=args[2]


# read input: 
input=read.table(input_file,header=F)


# setnames: 
colnames(input)=c('num_pc','num_peer','FDR0d1','FDR0d05','FDR0d01','FDR0d001')


# plot the number of egenes: 
p=ggplot(input,aes(x=as.factor(num_peer),y=FDR0d05,group=as.factor(num_pc),color=as.factor(num_pc)))+geom_point()+geom_line()+theme(axis.text.x=element_text(angle=90,vjust=0.5))+xlab('Number of PEER factors')+ylab('Number of eSNPs')+scale_color_discrete(name='Num PCs')+ggtitle('FDR=0.05')
save_plot(figure_path,p)