# Venn Diagram on Shared vs HCASMC-specific peak
# Boxiang Liu
# 2017-12-21
library(data.table)
library(VennDiagram)
library(cowplot)

fig_dir='../figures/hcasmc_specific_open_chromatin2/plot_venn_diagram/'
out_dir='../processed_data/hcasmc_specific_open_chromatin2/plot_venn_diagram/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir)}
if (!dir.exists(out_dir)) {dir.create(out_dir)}

intersect=fread('../processed_data/hcasmc_specific_open_chromatin2/intersect/intersect.bed')
hcasmc_col=which(names(intersect)=='HCASMC')
idx=which(intersect[,hcasmc_col,with=FALSE]==1)
n_hcasmc_peaks=length(idx)
n_hcasmc_spec_peaks=sum(intersect[idx,4,with=FALSE]==1)
n_shared=n_hcasmc_peaks-n_hcasmc_spec_peaks

data = data.frame(
	type = c('shared','hcasmc_specific'),
	value = c(n_shared, n_hcasmc_spec_peaks))

p=ggplot(data,aes(x='',y=value,fill=type,label=value))+
	geom_bar(width = 1, stat = "identity")+
	coord_polar('y',start=0)+
	annotate(geom='text',x='',y=n_hcasmc_peaks/2,label=paste(n_shared,sprintf('(%.02f%%)',n_shared/n_hcasmc_peaks*100)))+
	annotate(geom='text',x='',y=n_hcasmc_spec_peaks/2,label=paste(n_hcasmc_spec_peaks,sprintf('(%.02f%%)',n_hcasmc_spec_peaks/n_hcasmc_peaks*100)))+
	scale_fill_discrete(name='',labels=c('HCASMC','shared'))+
	theme(axis.text.x=element_blank(),
		axis.line=element_line(linetype='blank'),
		axis.title=element_text(size=0),
		axis.ticks=element_blank())
save_plot(sprintf('%s/pie_chart.pdf',fig_dir),p,base_height=4,base_width=4)
saveRDS(list(p=p),sprintf('%s/fig.rda',out_dir))

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