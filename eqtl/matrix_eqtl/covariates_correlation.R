library(data.table)
library(gplots)

# Variables:
in_fn='../processed_data/eqtl/matrix_eqtl/combine_covariates/covariates.tsv'
fig_dir='../figures/eqtl/matrix_eqtl/covariates_correlation/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

# Read data: 
covariates=read.table(in_fn,header=T,row.names=1,check.names=F)

cov_cor=cor(t(covariates))

pdf(sprintf('%s/covariates_correlation.pdf',fig_dir))
heatmap.2(cov_cor,Rowv=FALSE,Colv=FALSE,trace='none',margin=c(7,7))
dev.off()