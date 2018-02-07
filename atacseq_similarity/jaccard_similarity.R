library(stringr)
library(data.table)
library(cowplot)
library(ggrepel)
library(foreach)
library(doMC)
registerDoMC(20)

adult_dir='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_released_adult_tissue_group/'
fetal_dir='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_released_fetal_tissue_group/'
out_dir='../processed_data/atacseq_similarity/jaccard_similarity/'
fig_dir='../figures/atacseq_similarity/jaccard_similarity/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

sort_bed=function(in_dir,out_dir,pattern='bed'){
	if (!dir.exists(out_dir)){dir.create(out_dir,recursive=TRUE)}
	fn=list.files(in_dir,pattern=pattern,full.names=TRUE)
	foreach(i=fn)%dopar%{
		command=sprintf('sort -k1,1 -k2,2g "%s" > "%s/%s"',i,out_dir,basename(i))
		print(command)
		system(command)
	}
}

calc_jaccard=function(in_dir){
	fn=list.files(in_dir,pattern='bed',full.names=TRUE)
	hcasmc_fn=fn[str_detect(fn,'HCASMC')]

	jaccard=foreach(i = seq(length(fn)),.combine='rbind')%dopar%{
		sample=str_extract(fn[i],'(?<=sorted//)(.+?)(?=.merged)')
		print(sprintf('INFO - %s', sample))

		cmd=sprintf('bedtools jaccard -a "%s" -b "%s"',hcasmc_fn,fn[i])
		res=system(cmd,intern=TRUE)
		dt=data.table(t(as.numeric(str_split_fixed(res[2],'\t',4))))
		dt$sample=sample

		cmd=sprintf('bedtools jaccard -a "%s" -b "%s"',fn[i],fn[i])
		res=system(cmd,intern=TRUE)
		n_bases=str_split_fixed(res[2],'\t',4)[1]
		dt$n_bases=as.integer(n_bases)

		return(dt)
	}

	setnames(jaccard,c('intersection','union','jaccard','n_intersections','sample','n_bases'))
	setorder(jaccard,-jaccard)
	jaccard[,sample:=factor(sample,level=sample)]
	jaccard=jaccard[sample!='HCASMC',]
	jaccard$rank=seq(nrow(jaccard))

	return(jaccard)
}

add_gtex_name=function(jaccard){
	metadata_fn='../data/encode/dnase_seq/metadata.tsv'
	metadata=fread(metadata_fn)
	metadata=unique(metadata[`Biosample type`%in%c('tissue','primary cell')&`File Status`=='released'&`Biosample life stage`%in%c('fetal','adult'),list(gtex=`GTEx tissue`,sample=`Tissue group`)])
	jaccard=merge(jaccard,metadata,by='sample')
	return(jaccard)
}


get_color_map=function(){
	color=fread('shared/tissue_color.txt')
	color[,tissue_color_hex:=max(tissue_color_hex),by=tissue]
	color_map=color$tissue_color_hex
	names(color_map)=color$tissue_site_detail
	return(color_map)
}

# Adult samples:
sort_bed(adult_dir,out_dir=sprintf('%s/sorted/',adult_dir))
jaccard_adult=calc_jaccard(sprintf('%s/sorted/',adult_dir))
jaccard_adult=add_gtex_name(jaccard_adult)
jaccard_adult$life_stage='Adult'

# Fetal samples:
sort_bed(fetal_dir,out_dir=sprintf('%s/sorted/',fetal_dir))
jaccard_fetal=calc_jaccard(sprintf('%s/sorted/',fetal_dir))
jaccard_fetal=add_gtex_name(jaccard_fetal)
jaccard_fetal$life_stage='Fetal'
jaccard_fetal[,sample:=paste(sample,'(fetal)')]

# Merge adult and fetal:
jaccard=rbind(jaccard_adult,jaccard_fetal)
setorder(jaccard,-jaccard)
jaccard[,sample:=factor(sample,levels=sample)]
jaccard[,rank:=rank(-jaccard)]

# Get color map:
color_map=get_color_map()

# Make plot:
show_n=5
jaccard[,label:=ifelse( (rank<=show_n) | (rank>(length(jaccard)-show_n)),as.character(sample),'')]

top5_label=jaccard[rank<=show_n,as.character(sample)]
bottom5_label=jaccard[rank>(length(jaccard)-show_n),as.character(sample)]
bottom5_label=str_replace(bottom5_label,'inferior parietal ','')

p1=ggplot(jaccard,aes(sample,jaccard,label=label,color=gtex,shape=life_stage))+
	geom_point(size=2)+scale_y_log10(breaks=c(0,0.1,0.2,0.3),labels=c(0,0.1,0.2,0.3))+
	theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
	xlab(sprintf('ENCODE Tissue/Cell Type (n=%i)',nrow(jaccard)))+
	ylab('Epigenomic Similarity\n(Jaccard Index)')+
	annotate('text',x=10,y=0.05,label=paste(c('Top 5',top5_label),collapse='\n'))+
	annotate('text',x=49.5,y=0.13,label=paste(c('Bottom 5',bottom5_label),collapse='\n'))+
	scale_color_manual(values=color_map,guide='none')+
	scale_shape_discrete(name='')+
	theme(legend.position=c(0.8, 0.2))


# Save plot and data:
save_plot(sprintf('%s/jaccard_similarity.pdf',fig_dir),p1,base_aspect_ratio=2,base_height=3)
saveRDS(list(jaccard=jaccard,color_map=color_map,p1=p1),sprintf('%s/jaccard.rds',out_dir))


