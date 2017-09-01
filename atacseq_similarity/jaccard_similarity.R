library(stringr)
library(data.table)
library(cowplot)
library(ggrepel)
library(foreach)
library(doMC)
registerDoMC(20)

in_dir='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_filt/sorted/'
fig_dir='../figures/atacseq_similarity/jaccard_similarity'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

fn=list.files(in_dir,pattern='bed',full.names=TRUE)
hcasmc_fn=fn[str_detect(fn,'HCASMC')]



container=foreach(i = seq(length(fn)))%dopar%{

	sample=str_extract(fn[i],'(?<=sorted//)(.+?)(?=.merged)')

	print(sprintf('INFO - %s', sample))

	cmd=sprintf('bedtools jaccard -a "%s" -b "%s"',hcasmc_fn,fn[i])
	print(cmd)
	res=system(cmd,intern=TRUE)

	dt=data.table(t(as.numeric(str_split_fixed(res[2],'\t',4))))
	dt$sample=sample

	cmd=sprintf('bedtools jaccard -a "%s" -b "%s"',fn[i],fn[i])
	print(cmd)
	res=system(cmd,intern=TRUE)

	n_bases=str_split_fixed(res[2],'\t',4)[1]
	dt$n_bases=as.integer(n_bases)

	dt
}

jaccard=Reduce(rbind,container)

setnames(jaccard,c('intersection','union','jaccard','n_intersections','sample','n_bases'))
setorder(jaccard,-jaccard)
jaccard[,sample:=factor(sample,level=sample)]
jaccard=jaccard[sample!='HCASMC',]

jaccard$rank=seq(nrow(jaccard))
show_n=15
jaccard[,label:=ifelse( (rank<=show_n) | (rank>(length(jaccard)-show_n)),as.character(sample),'')]

p1=ggplot(jaccard,aes(sample,jaccard,label=label))+
	geom_point()+
	theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
	xlab(sprintf('ENCODE Sample (n=%i)',nrow(jaccard)))+
	ylab('Epigenomic Similarity (Jaccard Index)')+
	geom_text_repel(force=3,box.padding=unit(0.2, "lines"),nudge_x=1)

jaccard$sample
color=fread('shared/tissue_color.txt')

metadata_fn='../data/encode/dnase_seq/metadata.tsv'
metadata=fread(metadata_fn)
metadata=metadata[`Biosample type`%in%c('primary cell','tissue') & `Audit ERROR`=='']

jaccard=merge(jaccard,unique(metadata[,list(`Biosample term name`,`GTEx sample name`)]),by.x='sample',by.y="Biosample term name")
jaccard$GTEx
save_plot(sprintf('%s/jaccard_similarity.pdf',fig_dir),p1)
