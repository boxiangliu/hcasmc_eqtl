library(data.table)
library(foreach)
source('table1/utils.R')

sqtl_fn = '../processed_data/sqtl/fastQTL/adjust_pvalue/top_intron.txt'
out_dir = '../processed_data/table1/summarize_sqtl_discovery/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

read_sqtl_result = function(fn){
	sqtl = fread(fn,header=TRUE)
	return(sqtl)
}


sqtl = read_sqtl_result(sqtl_fn)
gene_annotation = read_gene_annotation(gene_annotation_fn)
sqtl = classify_genes(sqtl, gene_annotation)

condensed_sqtl_summary = foreach(FDR = c(1,0.05,0.01,0.001), .combine='rbind')%do%{
	significant_sqtl = sqtl[fdr<=FDR]
	sqtl_class_summary = summarize_genes(significant_sqtl)
	condensed_sqtl_summary = condense_summary(sqtl_class_summary, show = c('protein_coding','lincRNA','pseudogene'))
	condensed_sqtl_summary$fdr = FDR
	return(condensed_sqtl_summary)
}


reshaped_summary = reshape_summary(condensed_sqtl_summary)
setnames(reshaped_summary,'1','tested')
for(fdr in c('0.05','0.01','0.001')){
	setnames(reshaped_summary,fdr,'fdr')
	reshaped_summary$pct = calculate_percentage(reshaped_summary$fdr,reshaped_summary$tested)
	setnames(reshaped_summary,'fdr',fdr)
	setnames(reshaped_summary,'pct',sprintf('%s-pct',fdr))
}

out_path = sprintf('%s/sqtl_discovery_summary.txt',out_dir)
fwrite(reshaped_summary,out_path,sep='\t')
