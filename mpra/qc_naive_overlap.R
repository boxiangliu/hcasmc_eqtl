library(data.table)
library(dplyr)
library(cowplot)
args=commandArgs(T)
eqtl_file=args[1]
fig_dir=args[2]
eqtl_file='../processed_data/mpra/naive_overlap_eQTL/gwas_eQTL_naive_overlap.txt'
fig_dir='../figures/mpra/naive_overlap'


eqtl=fread(eqtl_file)
if (!dir.exists(fig_dir)) {dir.create(fig_dir)}


# plot rSNP r2:
p=ggplot(eqtl,aes(x=rsq_rsnp))+geom_density()+
geom_vline(xintercept=0.8,color='red',linetype=2)+
annotate(geom='text',x=0.6,y=30,label='ATP5G1,DHX36,TCF21,etc')+
geom_segment(aes(x=0.6,xend=0.88,y=28,yend=9),color='blue')+
xlab('Correlation b/t prior and posterior test SNP')
save_plot(paste0(fig_dir,'/rsq_rsnp.pdf'),p)


# plot rank:
cutoff=500
p1=ggplot(eqtl,aes(x=rank))+stat_bin(aes(y=cumsum(..count..)),geom='line',binwidth=2)+
geom_vline(xintercept=cutoff,color='red',linetype=2)+
annotate(geom='text',x=1000,y=2000,label=paste0(formatC(sum((eqtl$rank<cutoff)/nrow(eqtl))*100,digits=3),'% < ',cutoff))+
ylab('Cumulative density')+xlab('Rank')
save_plot(paste0(fig_dir,'/rank.pdf'),p1)


# plot effect size: 
p2=ggplot(eqtl[rsq_rsnp>0.8,],aes(x=pi))+geom_histogram()+
scale_x_continuous(breaks=seq(0,1,0.1))+
xlab('Allelic ratio')+ylab('Count')
save_plot(paste0(fig_dir,'/allelic_ratio.pdf'),p2)


# Plot number of variants vs p-value cutoff: 
cuts=c(1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8)
num_var=c()
for (cut in cuts){num_var=c(num_var,eqtl[pval<cut,rsid]%>%unique()%>%length())}
p3=ggplot(data.frame(cbind(cuts,num_var)),aes(-log10(cuts),num_var))+geom_point()+geom_text(label=num_var,nudge_y=50)+xlab('P-value cutoff')+ylab('Number of variants')
p4=ggplot(eqtl,aes(pval))+stat_ecdf()+xlab('P-value')+ylab('Empirical CDF')
pdf(paste0(fig_dir,'/pval_cut_vs_num_var.pdf'))
print(p3)
print(p4)
dev.off()