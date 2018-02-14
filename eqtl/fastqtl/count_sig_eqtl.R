# Count the number of significant eQTLs from fastQTL 
# permutation pass. 
# Boxiang Liu
# 2018-02-14

library(data.table)
library(foreach)
in_fn = '../processed_data/eqtl/fastqtl/output/permutation/all.txt.gz'
out_dir = '../processed_data/eqtl/fastqtl/output/permutation/'
if (!dir.exists(out_dir)) {dir.create(out_dir)}

read_eqtl = function(in_fn){
	eqtl = fread(sprintf('gunzip -c %s',in_fn),select=c(1,6,7,16),col.names=c('fid','sid','dist','pval'))
	return(eqtl)
}

adjust_pvalue = function(pval){
	padj = p.adjust(pval,method='fdr')
	return(padj)
}


count_sig_eqtl = function(padj,levels = c(0.05,0.01,0.001)){
	sig_eqtl = foreach (level = levels,.combine='rbind')%do%{
		n_sig = sum(padj < level, na.rm = TRUE)
		data.table(
			level = level,
			sig = n_sig
			)
	}
	return(sig_eqtl)
}

save_result = function(output, out_fn){
	fwrite(
		output,
		out_fn,
		sep = '\t'
		)
}
eqtl = read_eqtl(in_fn)
eqtl$padj = adjust_pvalue(eqtl$pval)
sig_eqtl = count_sig_eqtl(eqtl$padj)
save_result(sig_eqtl,sprintf('%s/sig_eqtl.tsv',out_dir))