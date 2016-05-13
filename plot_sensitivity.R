#!/usr/bin/env Rscript 
# bosh liu
# 2016/05/13
# durga
# plot sensitivity vs other metrics

# library:
library(cowplot)


# paths:
input_file="../data/joint2/recalibrate_INDEL.tranches"
figure_path='../figures/plot_sensitivity/'

# read input:
input=fread(input_file,header=T)


# remove sensitivity = 100%:
input=input[1:11,]


# plot sensitiviy vs numKnown, numNovel, minVQSLod:
p1=ggplot(input,aes(x=truthSensitivity,y=numKnown))+geom_point()+scale_x_continuous(breaks=round(input[,truthSensitivity],3))
p2=ggplot(input,aes(x=truthSensitivity,y=numNovel))+geom_point()+scale_x_continuous(breaks=round(input[,truthSensitivity],3))
p3=ggplot(input,aes(x=truthSensitivity,y=minVQSLod))+geom_point()+scale_x_continuous(breaks=round(input[,truthSensitivity],3))
p=plot_grid(p1,p2,p3,labels=c('A','B','C'),align='v',ncol=1)
save_plot(paste0(figure_path,'sensitivity.pdf'),p,base_height=12,base_width=6)