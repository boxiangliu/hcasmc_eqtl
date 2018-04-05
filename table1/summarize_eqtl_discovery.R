library(data.table)
library(foreach)
source('table1/utils.R')

tested_genes_fn = '../processed_data/rasqual/output_merged/expressed_genes.txt'
eqtl_genes_prefix = '../processed_data/rasqual/output_merged/treeQTL/eGenes_level'
out_dir = '../processed_data/table1/summarize_eqtl_discovery/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

read_tested_genes = function(fn){
	tested_genes = fread(fn,header=FALSE)
	setnames(tested_genes,'gene_id')
	return(tested_genes)
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




# Main: 
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

