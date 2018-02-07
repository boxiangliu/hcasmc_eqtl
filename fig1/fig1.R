library(data.table)
library(ggrepel)
library(cowplot)
library(ggplot2)
library(grid)
library(gridExtra)
library(VennDiagram)
library(stringr)
library(Gviz)
library(rtracklayer)
library(gtable)

out_dir='../figures/fig1/'
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# Fig. 1A:
load('../processed_data/mds/mds/mds.Rdata')
fig1a=p3

# Fig. 1B:
temp=readRDS('../processed_data/atacseq_similarity/jaccard_similarity/jaccard.rds')
fig1b=temp[['p1']]


# Fig. 1C:
temp=readRDS('../processed_data/intersect_NOISeq_DESeq2/fig.rda')
fig1c=temp[['p1b']]+theme(axis.text.x=element_text(angle=0,hjust=0.5))+scale_x_discrete(labels=c('CA','Fb'))



# Fig. 1D:
temp=readRDS('../processed_data/hcasmc_specific_open_chromatin2/plot_venn_diagram/fig.rda')
fig1d=temp[['p']]


# Fig. 1:
pdf(sprintf('%s/fig1.pdf',out_dir),height=7.5,width=15)
grid.arrange(arrangeGrob(fig1a,
	arrangeGrob(fig1b,
		arrangeGrob(fig1c,fig1d,ncol=2),
		arrangeGrob(rectGrob()),
		ncol=1),ncol=2))
dev.off()

setEPS()
postscript(sprintf('%s/fig1.eps',out_dir),height=7.5,width=15)
grid.arrange(arrangeGrob(fig1a,
	arrangeGrob(fig1b,
		arrangeGrob(fig1c,fig1d,ncol=2),
		arrangeGrob(rectGrob()),
		ncol=1),ncol=2))
dev.off()