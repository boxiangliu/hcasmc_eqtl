#!/usr/bin/env Rscript
# bosh liu
# durga
# 16/06/26
# plot verifyBamID "FREEMIX" column

# library:
library(cowplot)
library(data.table)

# paths: 
input_file='../processed_data/160526/detect_WGS_contamination/verifyBAMID.combined.tsv'
figure_path='../figures/160526/detect_WGS_contamination/'

# read input:
input=fread(input_file,header=F)


# setnames:
setnames(input, c('sample','contamination'))


# plot:
p=ggplot(input,aes(x=sample,y=contamination))+geom_point()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+geom_text(aes(label=ifelse(contamination>0.01,contamination,"")),angle=90,nudge_y=0.03)+background_grid(major = "x")+ylab('Contamination Proportion')+xlab('Sample')
save_plot(paste0(figure_path,'contamination_proportion.pdf'),p,base_width=8,base_height=8)