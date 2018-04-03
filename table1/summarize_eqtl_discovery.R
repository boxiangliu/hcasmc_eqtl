library(data.table)
library(foreach)

tested_genes_fn = '../processed_data/rasqual/output_merged/expressed_genes.txt'
gene_annotation_fn = '../data/gtex/gencode.v19.genes.v6p.hg19.bed'
eqtl_genes_prefix = '../processed_data/rasqual/output_merged/treeQTL/eGenes_level'
out_dir = '../processed_data/table1/summarize_eqtl_discovery/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

read_tested_genes = function(fn){
	tested_genes = fread(fn,header=FALSE)
	setnames(tested_genes,'gene_id')
	return(tested_genes)
}

read_gene_annotation = function(fn){
	gene_annotation = fread(fn)
	setnames(gene_annotation,c('chr','start','end','strand','gene_id','gene_name','type'))
	return(gene_annotation)
}

classify_genes = function(genes, gene_annotation){
	classified_genes = merge(genes, gene_annotation[,list(gene_id,type)], by = 'gene_id')
	return(classified_genes)
}

summarize_genes = function(classified_genes){
	class_summary = table(classified_genes$type)
	class_summary = data.table(class_summary)
	setnames(class_summary,c('type','num'))
	setorder(class_summary, -num)
	return(class_summary)
}

condense_summary = function(class_summary,show = c('protein_coding','lincRNA','pseudogene')){
	condensed_type = c()
	for (i in 1:nrow(class_summary)){
		type = class_summary$type[i]
		if (type %in% show){
			condensed_type = c(condensed_type,type)
		} else {
			condensed_type = c(condensed_type,'other')
		}
	}
	class_summary$condensed_type = condensed_type
	condensed_summary = class_summary[,list(num=sum(num)), by = 'condensed_type']
	return(condensed_summary)
}

read_eqtl_genes = function(fn){
	eqtl_genes = fread(fn)
	setnames(eqtl_genes,'family','gene_id') 
	return(eqtl_genes)
}


concat_eqtl_and_tested_gene_summaries = function(eqtl_summary,tested_gene_summary){
	tested_gene_summary$fdr = 'tested'
	concat_summary = rbind(eqtl_summary,tested_gene_summary)
	return(concat_summary)
}

reshape_summary = function(summary){
	reshaped_summary = dcast(summary,condensed_type~fdr,value.var='num')
	return(reshaped_summary)
}

calculate_percentage = function(x,y){
	percentage = as.numeric(sprintf('%0.2f',x/y*100))
	return(percentage)
}

tested_genes = read_tested_genes(tested_genes_fn)
gene_annotation = read_gene_annotation(gene_annotation_fn)
classified_genes = classify_genes(tested_genes, gene_annotation)
class_summary = summarize_genes(classified_genes)
condensed_summary = condense_summary(class_summary, show = c('protein_coding','lincRNA','pseudogene'))

condensed_eqtl_summary = foreach(fdr = c(0.05,0.01,0.001), .combine = 'rbind')%do%{
	eqtl_genes_fn = sprintf('%s%s.tsv',eqtl_genes_prefix,fdr)
	eqtl_genes = read_eqtl_genes(eqtl_genes_fn)
	classified_eqtl_genes = classify_genes(eqtl_genes, gene_annotation)
	eqtl_class_summary = summarize_genes(classified_eqtl_genes)
	condensed_eqtl_summary = condense_summary(eqtl_class_summary, show = c('protein_coding','lincRNA','pseudogene'))
	condensed_eqtl_summary$fdr = fdr
	return(condensed_eqtl_summary)
}

concat_summary = concat_eqtl_and_tested_gene_summaries(condensed_eqtl_summary,condensed_summary)
reshaped_summary = reshape_summary(concat_summary)
for(fdr in c('0.05','0.01','0.001')){
	setnames(reshaped_summary,fdr,'fdr')
	reshaped_summary$pct = calculate_percentage(reshaped_summary$fdr,reshaped_summary$tested)
	setnames(reshaped_summary,'fdr',fdr)
	setnames(reshaped_summary,'pct',sprintf('%s-pct',fdr))
}

out_path = sprintf('%s/eqtl_discovery_summary.txt',out_dir)
fwrite(reshaped_summary,out_path,sep='\t')