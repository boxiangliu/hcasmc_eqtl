library(data.table)
library(cowplot)

fig_dir='../figures/sqtl/optimal_covariate/compare_covariate/'
if (!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}

sig=fread('../processed_data/sqtl/optimal_covariate/test_covariate/sig.txt',col.names=c('FDR','sig','geno_pc','splice_pc'))

p=ggplot(sig[FDR%in%c(0.01,0.05)],aes(splice_pc,sig,color=as.character(geno_pc)))+
geom_point()+geom_line()+facet_grid(~FDR)+
scale_color_discrete(name='Genotype PC')+
xlab('Splice PC')+ylab('sQTL')

pdf(sprintf('%s/sig.pdf',fig_dir))
p
dev.off()