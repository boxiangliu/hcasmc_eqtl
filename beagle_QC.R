#!/urs/bin/env Rscript
# bosh liu
# 2016/05/15
# durga
# compare beagle Allelic R-Squared and Dosage R-Squared

# library:
library(cowplot)
library(data.table)
library(dplyr)

# paths:
input_file='../processed_data/160515_beagle_QC/recalibrated_biallelic_SNP.r2.tsv'
figure_path='../figures/160515_beagle_QC/'

# read input: 
input=fread(input_file,header=T)


# calculate the mean DR2 for each chromosome:
input1=input %>% group_by(CHROM) %>% summarize(DR2_mean=mean(DR2))
input1$CHROM=factor(input1$CHROM,levels=paste0('chr',seq(1,22)))


# plot R2 vs chrom:
p1=ggplot(input1,aes(x=CHROM,y=DR2_mean,group=1))+geom_line()+geom_point()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+xlab('Chromosome')+ylab('Dosage R2')
save_plot(paste0(figure_path,'dosageR2_by_chrom.pdf'),p1)


# bin AF:
breaks=seq(0,1,0.1)
input[,AF_bin:=cut(AF,breaks=breaks,include.lowest=T)]


# calculate mean dosage R2 by AF bins: 
input2=input%>%group_by(AF_bin)%>%summarize(DR2_mean=mean(DR2),n=n())
input2$AF_bin=factor(input2$AF_bin,levels=c('[0,0.1]','(0.1,0.2]','(0.2,0.3]','(0.3,0.4]','(0.4,0.5]','(0.5,0.6]','(0.6,0.7]','(0.7,0.8]','(0.8,0.9]','(0.9,1]'))


# make histogram of Allelic R2 binned by AF:
p2=ggplot(input2,aes(x=AF_bin,y=DR2_mean,group=1))+geom_line()+geom_point()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+xlab('AF')+ylab('Dosage R2')
p3=ggplot(input2,aes(x=AF_bin,y=n,group=1))+geom_line()+geom_point()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+xlab('AF')+ylab('Number of sites')
p4=plot_grid(p2,p3,labels=c('A','B'))
save_plot(paste0(figure_path,'dosageR2_by_AFbin.pdf'),p4,base_width=8)


# bin AF:
breaks=seq(0,1,0.002)
input[,AF_bin:=cut(AF,breaks=breaks,include.lowest=T)]


# calculate mean dosage R2 by AF bins: 
input2=input%>%group_by(AF_bin)%>%summarize(DR2_mean=mean(DR2),n=n())
input2[,AF_lb:=str_split_fixed(AF_bin,"\\(|\\[|,",3)[,2]%>%as.numeric()]


# make histogram of Allelic R2 binned by AF:
p2=ggplot(input2,aes(x=AF_lb,y=DR2_mean,group=1))+geom_line()+geom_point()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+xlab('AF')+ylab('Dosage R2')
p3=ggplot(input2,aes(x=AF_lb,y=n,group=1))+geom_line()+geom_point()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+xlab('AF')+ylab('Number of sites')
p4=plot_grid(p2,p3,labels=c('A','B'))
save_plot(paste0(figure_path,'dosageR2_by_AFbin_size_002.pdf'),p4,base_width=8)


# plot histogram of R2:
pdf(paste0(figure_path,'R2_histogram.pdf'))
histogram=hist(input$DR2,breaks=seq(0,1,0.1),ylim=c(0,2e7),main='Dosage R2',xlab='Dosage R2')
with(histogram,text(mids,counts+2e6,paste0(round(10*density,3),'%'),srt=90))
dev.off()
