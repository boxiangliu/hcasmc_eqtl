#!/usr/bin/env Rscript
# bosh liu
# 2016/05/11
# durga
# plot allele frequency of pruned marker set



input_file='../processed_data/160511_calc_freq/indep_pairwise_50_5_0.2.frq'
input=fread(input_file)


# make histogram of allele frequency distribution:
pdf('../figures/160511_calc_freq/maf_dist.pdf')
hist(input$MAF,main='MAF of pruned set',xlab='MAF')
abline(v=0.05,col='red',lty=2)
dev.off()