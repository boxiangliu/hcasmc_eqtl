#!/usr/bin/env Rscript
# bosh liu
# 2016/05/19
# durga
# calculate RNAseq sample-sample correlation using rpkm values


# libraries:
library(data.table)
library(gplots)
library(cowplot)


# paths:
input_file='../processed_data/160519_rpkm/combined.rpkm'
figure_path='../figures/160519_rpkm/'


# read rpkm:
rpkm=fread(input_file,header=T)


# remove the first two columns: 
rpkm2=rpkm
rpkm2[,Name:=NULL]
rpkm2[,Description:=NULL]


# calculate correlation:
cor_mat=cor(rpkm2)
pdf(paste0(figure_path,'rpkm_correlation_heatmap.pdf'))
heatmap.2(cor_mat,Rowv=F,Colv=F,trace='none',dendrogram='none')
dev.off()




# calculate mean correlation of each sample with other samples:
diag(cor_mat)=NaN
mean_cor=rowMeans(cor_mat,na.rm=T)


# calculate overall mean correlation: 
overall_mean_cor=mean(cor_mat,na.rm=T)


# calculate D-statistic: 
mad=median(abs(mean_cor-overall_mean_cor))
D=(mean_cor-overall_mean_cor)/mad


# plot:
to_plot=data.frame(sample=names(D),D=D)
p=ggplot(to_plot,aes(x=sample,y=D))+geom_point()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+background_grid(major = "x", minor = "x")
save_plot(paste0(figure_path,'D_statistic.pdf'),p,base_width=8,base_height=8)