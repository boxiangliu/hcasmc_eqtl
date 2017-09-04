library(data.table)
library(ggrepel)
library(cowplot)
library(ggplot2)
library(grid)
library(gridExtra)

out_dir='../figures/fig1/'
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# Fig. 1A:
load('../processed_data/160603/mds.Rdata')
fig1a=p3


# Fig. 1B:
temp=readRDS('../processed_data/atacseq_similarity/jaccard_similarity/jaccard.rds')
jaccard=temp[[1]]
color_map=temp[[2]]


show_n=5
jaccard[,label:=ifelse( (rank<=show_n) | (rank>(length(jaccard)-show_n)),as.character(sample),'')]


fig1b=ggplot(jaccard,aes(sample,jaccard,label=label,color=gtex))+
	geom_point(size=1)+
	theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
	xlab(sprintf('ENCODE Sample (n=%i)',nrow(jaccard)))+
	ylab('Epigenomic Similarity\n(Jaccard Index)')+
	geom_text_repel(force=5)+
	scale_color_manual(values=color_map,guide='none')



# Fig. 1:
pdf(sprintf('%s/fig1.pdf',out_dir),height=7,width=12)
grid.arrange(fig1a,arrangeGrob(fig1b,fig1b,fig1b,ncol=1),ncol=2,widths=c(3,3))
dev.off()