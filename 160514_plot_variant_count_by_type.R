#!/usr/bin/env Rscript
# bosh liu
# durga
# 2016/05/14
# plot the number of variants by type

# library: 
library(data.table)
library(dplyr)
library(cowplot)


# paths:
input_file='../processed_data/160514_variant_count_by_type/count.txt'
figure_path='../figures/160514_variant_count_by_type/'

# read input:
input=fread(input_file,header=F)


# subset to useful columns:
input[,V1:=NULL]
input[,V2:=NULL]


# setnames:
setnames(input, c('type','count'))


# plot: 
input[,type:=str_replace(type,'number of ','') %>% str_replace(':','')]
input=rbind(input,data.table(type='biallelic SNP',count=input[type=='SNPs',count]-input[type=='multiallelic SNP sites',count]))
input2=input[c(2,3,9,5)]
input2$type=factor(input2$type,levels=c('records','SNPs','biallelic SNP','indels'))
p=ggplot(input2,aes(x=type,y=count,label=count))+geom_bar(stat='identity')+geom_text(aes(y=count+5e5))
save_plot(paste0(figure_path,'variant_count_by_type.pdf'),p)

