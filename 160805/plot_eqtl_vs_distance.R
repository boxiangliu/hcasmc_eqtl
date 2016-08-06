# plot the number of significant eqtl vs distance to intron start site. 


# command args:
args=commandArgs(T,T)
eqtl_file=args$input
num_eqtl_vs_dist_fig=args$num_eqtl_vs_dist
pval_vs_dist_fig=args$pval_vs_dist

# eqtl_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/hcasmc.eqtl.pc4.peer8.padj.txt'
# num_eqtl_vs_dist_fig='/srv/persistent/bliu2/HCASMC_eQTL/figures/160805/num_sig_eqtl_vs_dist.pdf'
# pval_vs_dist_fig='/srv/persistent/bliu2/HCASMC_eQTL/figures/160805/eqtl_pval_vs_dist.pdf'

# read eqtl:
eqtl=fread(eqtl_file)


# subset to significant eqtls:
sig_eqtl=eqtl[qval<0.05,]


# plot number of sig. eqtls vs distance:
p=ggplot(sig_eqtl,aes(dist))+geom_histogram(binwidth=50)
save_plot(num_eqtl_vs_dist_fig,p)


# plot eqtl pvalue againt distance to TSS: 
eqtl[,log_pval:=-log10(pval)]
set.seed(42)
idx=sample(1:nrow(eqtl),1e5)
p1=ggplot(eqtl[idx,],aes(x=dist,y=log_pval,size=log_pval))+geom_point(alpha=0.3)+xlab('Distance to TSS')+ylab('-log10(p-value)')
save_plot(pval_vs_dist_fig,p1)