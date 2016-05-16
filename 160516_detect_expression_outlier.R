#!/usr/bin/env Rscript
# bosh liu
# durga
# 2016/05/16
# calculate sample-sample correlation and detect outliers

# library:
library(cowplot)


# paths: 
input_file='../processed_data/030_variance_stabilize/counts_size_corrected.rds'
figure_path='../figures/160516_detect_expression_outlier/'

# read size factor corrected counts:
input=readRDS(input_file)


# calculate correlation between all samples: 
correlation=cor(input)
diag(correlation)=NaN


# calculate mean correlation for each sample:
mean_cor=rowMeans(correlation,na.rm=T)


# calculate overall mean correlation: 
overall_mean_cor=mean(correlation,na.rm=T)


# calculate D-statistic: 
mad=median(abs(mean_cor-overall_mean_cor))


# calculate D statistic:
D=(mean_cor-overall_mean_cor)/mad


# plot D statistic:
to_plot=data.frame(sample=names(D),D=D)
p=ggplot(to_plot,aes(x=sample,y=D))+geom_point()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+background_grid(major = "x", minor = "x")
save_plot(paste0(figure_path,'D_statistic.pdf'),p,base_width=8,base_height=8)
