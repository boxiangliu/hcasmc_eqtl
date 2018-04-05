gene_annotation_fn = '../data/gtex/gencode.v19.genes.v6p.hg19.bed'

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


reshape_summary = function(summary){
	reshaped_summary = dcast(summary,condensed_type~fdr,value.var='num')
	return(reshaped_summary)
}

calculate_percentage = function(x,y){
	percentage = as.numeric(sprintf('%0.2f',x/y*100))
	return(percentage)
}