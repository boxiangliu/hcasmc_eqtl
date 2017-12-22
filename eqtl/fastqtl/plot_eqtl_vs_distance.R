# plot the number of significant eqtl vs distance to intron start site. 
library(R.utils)
library(data.table)
library(stringr)
library(cowplot)

# command args:
args=commandArgs(T,T)
eqtl_file=args$eqtl_file
pval_vs_dist_fig=args$pval_vs_dist

eqtl_file='../processed_data/eqtl/fastqtl/output/nominal/all.txt.gz'
pval_vs_dist_fig='../figures/eqtl/fastqtl/plot_eqtl_vs_distance/eqtl_pval_vs_dist.pdf'

# read eqtl:
if (str_detect(eqtl_file,'.gz')) {
	eqtl=fread(sprintf('zcat %s',eqtl_file),header=F)
} else {
	eqtl=fread(eqtl_file,header=F)
}

setnames(eqtl,c('fid','snp','distance','pval','beta','varbeta'))


# plot eqtl pvalue againt distance to TSS: 
eqtl[,log_pval:=-log10(pval)]
set.seed(42)
idx=sample(1:nrow(eqtl),1e5)
p1=ggplot(eqtl[idx,],aes(x=distance,y=log_pval,size=log_pval))+geom_point(alpha=0.3)+xlab('Distance to TSS')+ylab('-log10(p-value)')+scale_size(guide='none')
save_plot(pval_vs_dist_fig,p1)