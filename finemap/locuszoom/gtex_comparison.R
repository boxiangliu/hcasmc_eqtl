library(data.table)
library(stringr)
library(cowplot)
library(ggrepel)

gtex_fm_fn='../processed_data/finemap/finemap_mike/sorted_gtex_finemap_hcasmc_results.txt'
hcasmc_fm_fn='../processed_data/finemap/finemap_mike/UKBB_GWAS1KG_EXOME_CAD_SOFT_META_PublicRelease_300517_txt_gz_finemap_clpp_status.txt'
fig_dir='../figures/finemap/locuszoom/gtex_comparison/'
if (!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}

tissue_id2name=function(tissue_id){
	color=fread('shared/tissue_color.txt')
	color[tissue_site_detail_id=='hcasmc.rpkm',tissue_site_detail_id:='hcasmc']
	map=color$tissue_site_detail
	names(map)=color$tissue_site_detail_id
	return(map[tissue_id])
}

read_finemap=function(fm_fn,threshold=5e-5,subset=NULL){

	tmp=fread(fm_fn,col.names=c('snp','eqtl_file','gwas_file','gene_name','conditional_level','n_tested_snps','clpp_score','gwas_logp','eqtl_logp'))
	tmp[,chrom:=str_split_fixed(snp,'_',2)[,1]]
	tmp[,pos:=as.integer(str_split_fixed(snp,'_',2)[,2])]
	fm=tmp[grepl('eqtl',eqtl_file),list(chrom,pos,tissue=eqtl_file,gene_name,y=clpp_score,gwas_logp,eqtl_logp,method='eCAVIAR')]

	set.seed(42)
	fm[,rank:=rank(-y,ties.method='random'),by=c('gene_name','tissue')]
	fm=fm[rank==1]
	fm$rank=NULL
	fm=fm[gwas_logp> -log10(threshold)|eqtl_logp> -log10(threshold)]
	fm[,tissue:=gsub('_Analysis_cis_eqtl_gz|_eqtls_txt_gz','',tissue)]
	fm$tissue=tissue_id2name(fm$tissue)
	if (!is.null(subset)) {
		fm=fm[gene_name%in%subset]
	}
	return(fm)
}

gene_id2name=function(id,map){
	return(map[id])
}

abbreviate_tissue_name=function(tissue_name){
	color=fread('shared/tissue_color.txt')
	color[tissue_site_detail_id=='hcasmc.rpkm',tissue_site_detail_id:='hcasmc']
	map=color$abbreviation
	names(map)=color$tissue_site_detail
	return(map[tissue_name])
}


get_color_map=function(){
	color=fread('shared/tissue_color.txt')
	color[,tissue_color_hex:=max(tissue_color_hex),by=tissue]
	color_map=color$tissue_color_hex
	names(color_map)=color$tissue_site_detail
	return(color_map)
}


plot_clpp=function(data,color_map,top=length(unique(data$tissue)),label_top=length(unique(data$tissue))){
	setorder(data,-y)
	data=data[1:top]
	setorder(data,y)
	data[,tissue:=factor(tissue,tissue)]
	data[,abbreviation:=factor(abbreviation,abbreviation)]
	data[,label:=ifelse(rank(-y)<=label_top,abbreviation,'')]
	p=ggplot(data,aes(y,abbreviation,color=tissue,label=label))+
		geom_point()+
		scale_color_manual(values=color_map,guide='none')+
		xlab(sprintf('%s CLPP',unique(data$gene_name)))+ylab('GTEx tissues')+ 
		theme(axis.text.y=element_text(color=ifelse(data$abbreviation=='HCASMC','purple','black')))
		#+
		# theme(axis.text.y=element_blank(),
			# axis.ticks.y=element_blank())
	return(p)
}



gtex_fm=read_finemap(gtex_fm_fn,threshold=1)
gtex_fm$gene_name=gene_id2name(gtex_fm$gene_name,c(ENSG00000118526.6='TCF21',ENSG00000134853.7='PDGFRA',ENSG00000166949.11='SMAD3',ENSG00000182511.7='FES',ENSG00000213445.4='SIPA1'))
hcasmc_fm=read_finemap(hcasmc_fm_fn,threshold=1,subset=c('TCF21','SIPA1','SMAD3','FES','PDGFRA'))
fm=rbind(gtex_fm,hcasmc_fm)
fm$abbreviation=abbreviate_tissue_name(fm$tissue)
color_map=get_color_map()

pdf(sprintf('%s/gtex_comparison.pdf',fig_dir),height=6,width=6)
for (i in unique(fm$gene_name)){
	data=fm[gene_name==i]
	p=plot_clpp(data,color_map)
	print(p)
}
dev.off()

pdf(sprintf('%s/gtex_comparison.top20.pdf',fig_dir),height=4,width=5)
for (i in unique(fm$gene_name)){
	data=fm[gene_name==i]
	p1=plot_clpp(data,color_map,top=20)+geom_point(size=2)
	print(p1)
}
dev.off()