library(data.table)
library(cowplot)
library(gplots)

# Variables:
in_fn='../processed_data/eqtl/peer/factors.tsv'
fig_dir='../figures/eqtl/peer/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

# Read data:
factors=fread(in_fn,header=T)
setDF(factors)
rownames(factors)=factors$ID
factors$ID=NULL

fcor=cor(t(factors))

pdf(sprintf('%s/peer_factor_correlation.pdf',fig_dir))
heatmap.2(fcor,Rowv=FALSE,Colv=FALSE,trace='none')
dev.off()