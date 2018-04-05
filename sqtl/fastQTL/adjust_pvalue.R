library(data.table)
library(stringr)
library(leafcutter)

sqtl_fn = '../processed_data/sqtl/fastQTL/permutation/all.permutation.txt.gz'
gtf_fn = '/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6_v6p_annotations/gencode.v19.genes.v6p.patched_contigs.gtf.gz'
exon_file = '/srv/persistent/bliu2/tools/leafcutter/leafcutter/data/gencode19_exons.gene_id.txt.gz'
expressed_genes_fn = '../processed_data/rasqual/output_merged/expressed_genes.txt'
gene_annotation_fn = '../data/gtex/gencode.v19.genes.v6p.hg19.bed'
out_dir = '../processed_data/sqtl/fastQTL/adjust_pvalue/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

gtf_to_exon = function(gft_fn){
	gtf = read.table(gtf_fn,stringsAsFactors=FALSE,header=FALSE,sep='\t',quote="")[,c(1,3,4,5,7,9)]
	setDT(gtf)
	setnames(gtf,c('chr','feature','start','end','strand','annotation'))
	exon = gtf[feature == 'exon']
	exon[,chr := paste0('chr',chr)]
	exon[,gene_name := str_extract(annotation,'(?<=gene_id\\s\\")(ENSG.+?)(?=\\";)')]
	exon$feature = NULL
	exon$annotation = NULL
	return(exon)
}


read_expressed_genes = function(fn){
	expressed_genes = unlist(fread(fn,header=FALSE))
	return(expressed_genes)
}

subset_to_expressed_genes = function(sqtl,expressed_genes){
	expressed_gene_sqtl = sqtl[gene_id %in% expressed_genes]
	return(expressed_gene_sqtl)
}

read_fastqtl_permutation_mode = function(fn){
	fastqtl = read.table(sqtl_fn,header=FALSE,stringsAsFactors=FALSE)[,c(1,16)]
	setDT(fastqtl)
	setnames(fastqtl,c('intron','pval'))
	return(fastqtl)
}

extract_intron_cluster = function(intron){
	intron_cluster = str_split_fixed(intron,':',4)[,4]
	return(intron_cluster)
}

bonferroni_correction_by_intron_cluster = function(x){
	x = copy(x)
	x[,n := .N, by = 'cluster']
	x[,bonf := pval*n]
	x[,bonf := ifelse(bonf > 1, 1, bonf)]
	return(x$bonf)
}

select_top_intron_per_cluster = function(sqtl){
	sqtl[, bonf_rank := rank(bonf,ties.method = 'first'), by = 'cluster']
	top_intron = sqtl[bonf_rank == 1]
	top_intron$bonf_rank = NULL
	return(top_intron)
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

save_top_intron = function(top_intron,fn){
	top_intron$cluster = NULL
	fwrite(top_intron,fn,sep='\t')
}

exon = gtf_to_exon(gtf_fn)
write.table(exon,gzfile(exon_file),col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')

sqtl = read_fastqtl_permutation_mode(sqtl_fn)
sqtl$gene_id=assign_gene(exon_file,sqtl$intron)
expressed_genes = read_expressed_genes(expressed_genes_fn)
sqtl = subset_to_expressed_genes(sqtl,expressed_genes)

sqtl$cluster = extract_intron_cluster(sqtl$intron)
sqtl$bonf = bonferroni_correction_by_intron_cluster(sqtl[,list(cluster,pval)])

top_intron = select_top_intron_per_cluster(sqtl)
top_intron$fdr = p.adjust(top_intron$bonf, method='fdr')

out_fn = paste0(out_dir,'/top_intron.txt')
save_top_intron(top_intron,out_fn)
