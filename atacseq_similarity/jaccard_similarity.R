library(stringr)
library(data.table)
library(cowplot)
library(ggrepel)
library(foreach)
library(doMC)
registerDoMC(20)

in_dir='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks/'
in_dir_filt='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_filt/sorted/'
in_dir_released='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_released/'
out_dir='../processed_data/atacseq_similarity/jaccard_similarity/'
fig_dir='../figures/atacseq_similarity/jaccard_similarity'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

calc_jaccard=function(in_dir){

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

	jaccard
}

jaccard=calc_jaccard(in_dir_filt)


metadata_fn='../data/encode/dnase_seq/metadata.tsv'
metadata=fread(metadata_fn)
metadata=metadata[`Biosample type`%in%c('primary cell','tissue')&`Audit ERROR`=="",]

jaccard=merge(jaccard,unique(metadata[,list(`Biosample term name`,gtex=`GTEx tissue`)]),by.x='sample',by.y="Biosample term name")
stopifnot(nrow(jaccard)==102)

color=fread('shared/tissue_color.txt')
color[,tissue_color_hex:=max(tissue_color_hex),by=tissue]

jaccard[,gtex:=ifelse(gtex%in%color$tissue,gtex,'Missing')]

temp=unique(color[,list(tissue,tissue_color_hex)])
color_map=c(temp$tissue_color_hex,'#000000')
names(color_map)=c(temp$tissue,'Missing')
setorder(jaccard,-jaccard)
jaccard[,sample:=factor(sample,levels=sample)]

show_n=15
jaccard[,label:=ifelse( (rank<=show_n) | (rank>(length(jaccard)-show_n)),as.character(sample),'')]
saveRDS(list(jaccard=jaccard,color_map=color_map),sprintf('%s/jaccard.rds',out_dir))
# jaccard=readRDS(sprintf('%s/jaccard.rds',out_dir))[[1]]

p1=ggplot(jaccard,aes(sample,jaccard,label=label,color=gtex))+
	geom_point()+
	theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
	xlab(sprintf('ENCODE Sample (n=%i)',nrow(jaccard)))+
	ylab('Epigenomic Similarity\n(Jaccard Index)')+
	geom_text_repel(force=3,box.padding=unit(0.2, "lines"),nudge_x=1)+
	scale_color_manual(values=color_map,guide='none')


save_plot(sprintf('%s/jaccard_similarity.pdf',fig_dir),p1,base_aspect_ratio=3,base_height=3)



