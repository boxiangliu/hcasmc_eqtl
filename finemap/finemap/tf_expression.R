library(data.table)
library(cowplot)

fig_dir='../figures/finemap/finemap/tf_expression/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

# Read count data:
count_fn='../data/rnaseq2/read_count/rnaseqc/rnaseqc.hcasmc_eqtl.reads.gct'
count=fread(count_fn,header=TRUE)

# TCF21 - rs2327429:
tf=c('THAP1', 'CENPB', 'NFIX')
tf_expression=melt(count[Description%in%tf,2:ncol(count)],id.vars='Description',variable.name='sample',value.name='count')
p1=ggplot(tf_expression,aes(Description,count))+geom_boxplot(outlier.size=-1)+geom_jitter(width=0.1)+scale_y_log10()+xlab('Gene')+ylab('Count')
save_plot(sprintf('%s/tcf21_rs2327429.pdf',fig_dir),p1)