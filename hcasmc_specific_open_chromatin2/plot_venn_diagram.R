# Venn Diagram on Shared vs HCASMC-specific peak
# Boxiang Liu
# 2017-12-21
library(data.table)
library(VennDiagram)
fig_dir='../figures/hcasmc_specific_open_chromatin2/plot_venn_diagram/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir)}

intersect=fread('../processed_data/hcasmc_specific_open_chromatin2/intersect/intersect.bed')
hcasmc_col=which(names(intersect)=='HCASMC')
idx=which(intersect[,hcasmc_col,with=FALSE]==1)
n_hcasmc_peaks=length(idx)
n_hcasmc_spec_peaks=sum(intersect[idx,4,with=FALSE]==1)

pdf(sprintf('%s/venn_diagram.pdf',fig_dir),height=4,width=4)
fig1e=draw.pairwise.venn(
	area1=n_hcasmc_peaks, 
	area2=nrow(intersect),
	cross.area=n_hcasmc_peaks-n_hcasmc_spec_peaks,
	category=c("HCASMC", "ENCODE"),
	lty=rep("blank", 2), 
	fill=c("blue", "red"), 
	alpha = rep(0.5, 2), 
	cat.pos = c(0, 0), 
	cat.dist = rep(0.025, 2),
	cex=1.3,
	cat.cex=1.5)
dev.off()