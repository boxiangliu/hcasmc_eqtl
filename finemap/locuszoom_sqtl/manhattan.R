library(data.table)
library(cowplot)
library(ggrepel)
library(manhattan)
library(stringr)
library(leafcutter)

smr_fn='../processed_data/finemap/smr/UKBB_HCASMC_sqtl_smr_results5e-05.out.smr'
fm_fn='../processed_data/finemap/finemap_mike/UKBB_GWAS1KG_EXOME_CAD_SOFT_META_PublicRelease_300517_txt_gz_finemap_clpp_status.txt'
exon_file='/srv/persistent/bliu2/tools/leafcutter/leafcutter/data/gencode19_exons.txt.gz'
fig_dir='../figures/finemap/locuszoom_sqtl/manhattan/'
out_dir='../processed_data/finemap/locuszoom_sqtl/manhattan/'
if(!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}
if(!dir.exists(out_dir)){dir.create(out_dir,recursive=TRUE)}

read_smr=function(smr_fn){
	fread(smr_fn)[,list(chrom=ProbeChr,clu_name=paste(Gene,cluster,sep='_'),pos=Probe_bp,y=log10(p_SMR),gwas_logp=-log10(p_GWAS),sqtl_logp=-log10(p_eQTL),method='SMR')]
}

read_finemap=function(fm_fn,threshold=5e-5){
	if (grepl('clpp_status',fm_fn)) {
		tmp=fread(fm_fn,col.names=c('snp','sqtl_file','gwas_file','clu_name','conditional_level','n_tested_snps','clpp_score','gwas_logp','sqtl_logp'))
		tmp[,chrom:=str_split_fixed(snp,'_',2)[,1]]
		tmp[,pos:=as.integer(str_split_fixed(snp,'_',2)[,2])]
		fm=tmp[grepl('sqtl',sqtl_file),list(chrom,pos,clu_name,y=clpp_score,gwas_logp,sqtl_logp,method='eCAVIAR')]
	} else {
		fm=fread(fm_fn)[,list(chrom=chrom,pos,clu_name=gene,y=clpp_score,gwas_logp=gwas_log_pval,sqtl_logp=sqtl_log_pval,method='eCAVIAR')]
	}

	set.seed(42)
	fm[,rank:=rank(-y,ties.method='random'),by='clu_name']
	fm=fm[rank==1]
	fm$rank=NULL
	fm=fm[gwas_logp> -log10(threshold)|sqtl_logp> -log10(threshold)]
	return(fm)
}


assign_gene=function(exon_file,clusters){
	exons_table     = read.table(exon_file, header=T, stringsAsFactors = F)
	intron_meta     = get_intron_meta(clusters)
	exons_table$chr = add_chr(exons_table$chr)
	intron_meta$chr = add_chr(intron_meta$chr)
	clu_gene_map    = map_clusters_to_genes(intron_meta, exons_table)
	clu_map = data.frame(clusters,clu=gsub(':[0-9]+:[0-9]+','',clusters))
	clu_gene_map = merge(clu_gene_map,clu_map,by='clu')
	map = clu_gene_map$genes
	names(map) = clu_gene_map$clusters
	return(map[clusters])
}

subset_to_protein_coding_and_lncRNA=function(x){
	annotation=fread(
		input = '../data/gtex/gencode.v19.genes.v6p.hg19.bed',
		select = c(6,7),
		col.names = c('gene_name','type')
		)
	annotation = annotation[type %in% c('protein_coding','lincRNA')]
	x = x[gene_name %in% annotation$gene_name]
	return(x)
}

smr=read_smr(smr_fn)
smr$clusters=paste0('chr',gsub('clu:','clu_',gsub('_',':',smr$clu_name)))
smr$gene_name=assign_gene(exon_file,smr$clusters)
smr=subset_to_protein_coding_and_lncRNA(smr)


fm=read_finemap(fm_fn,threshold=1e-4)
fm$clusters=paste0('chr',gsub('clu:','clu_',gsub('_',':',fm$clu_name)))
fm$gene_name=assign_gene(exon_file,fm$clusters)
fm=subset_to_protein_coding_and_lncRNA(fm)


smr_threshold=0.05/2439
fm_threshold=0.05
data=rbind(smr,fm)
data[,method:=factor(method,level=c('SMR','eCAVIAR'))]
data[,label:=ifelse( (y<log10(smr_threshold)) | (y>fm_threshold),paste(gene_name,str_split_fixed(clusters,':clu_',2)[,1],sep='\n'),'')]
data[,chrom:=paste0('chr',chrom)]

dummy=data.table(method=c('SMR','eCAVIAR'),y=c(log10(smr_threshold),fm_threshold))

p=manhattan(data,build='hg19')+
	facet_grid(method~.,scale='free_y')+
	scale_y_continuous(labels=function(x){abs(x)})+
	geom_text_repel(aes(label=label),force=3)+
	geom_hline(data=dummy,aes(yintercept=y),color='red',linetype=2)+
	ylab(paste('-log10(P)                 CLPP'))

p=manhattan(data,build='hg19')+
	facet_grid(method~.,scale='free_y')+
	scale_y_continuous(labels=function(x){abs(x)})+
	geom_text(aes(label=label),hjust=1.1)+
	geom_hline(data=dummy,aes(yintercept=y),color='red',linetype=2)+
	ylab(paste('-log10(P)                 CLPP'))


saveRDS(list(data,dummy,p),sprintf('%s/manhattan.rds',out_dir))
save_plot(sprintf('%s/manhattan.pdf',fig_dir),p,base_width=8,base_height=4)


