# devtools::install_github('dviraran/xCell')
library(xCell)
library(GSEABase)
library(data.table)
library(foreach)
source('/srv/persistent/bliu2/rpe/scripts/utils/genome_annotation.R')
library(cowplot)

hcasmc_signature_fn = '../processed_data/gwas_gene_overlap/ldscore_regression/tissue_specific_gene/tissue_specific_gene_no_sm/HCASMC.3sd.txt'
coronary_artery_read_count_fn = '/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6p_read_counts/Artery_Coronary.reads.txt'
fig_dir = '../figures/rebuttal/reviewer1/cell_type_deconvolution/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

read_gene_list = function(fn){
	x = fread(fn)
	return(x$Description)
}

get_xCell_signatures = function(){
	signatures = foreach(i = 1:length(xCell.data$signature),.combine='c')%dopar%{
		gene_ids = geneIds(xCell.data$signature[[i]])
		name = setName(xCell.data$signature[[i]])
		x = list(gene_ids)
		names(x) = name
		return(x)
	}
	return(signatures)
}

read_coronary_artery_read_count = function(coronary_artery_read_count_fn){
	coronary_artery_read_count = fread(coronary_artery_read_count_fn) 
	coronary_artery_read_count = merge(gene_annotation[,list(gene_id,gene_name)],coronary_artery_read_count, by.x = 'gene_id', by.y = 'Gene', sort = FALSE)
	setDF(coronary_artery_read_count)
	rownames(coronary_artery_read_count) = coronary_artery_read_count$gene_name
	coronary_artery_read_count = coronary_artery_read_count[,-c(1,2)]
	return(coronary_artery_read_count)
}

make_plot_data = function(res){
	res.t = as.data.frame(t(res))
	res_wide = melt(res.t,variable.name='cell_type',value.name='proportion')
	setDT(res_wide)
	res_wide_median = res_wide[,list(median=median(proportion)),by='cell_type']
	setorder(res_wide_median,-median)
	res_wide = res_wide[cell_type%in%res_wide_median[1:30,cell_type],]
	res_wide = res_wide[!cell_type%in%c('MicroenvironmentScore','StromaScore'),]
	res_wide[,cell_type:=factor(cell_type,levels=rev(res_wide_median$cell_type))]
	return(res_wide)
}

plot_xCell = function(res_wide){
	p = ggplot(res_wide,aes(x=cell_type,y=proportion)) + 
		geom_boxplot(outlier.size=-1) + 
		xlab('') + 
		ylab('xCell percentages') + 
		coord_flip()
	return(p)
}

hcasmc_signature = read_gene_list(hcasmc_signature_fn)
signatures = get_xCell_signatures()
signatures[['HCASMC']] = hcasmc_signature

gene_annotation = read_gencode()[,list(gene_id,gene_name)]
gene_annotation = gene_annotation[!duplicated(gene_name),]
coronary_artery_read_count = read_coronary_artery_read_count(coronary_artery_read_count_fn)

res = xCellAnalysis(coronary_artery_read_count, signature=signatures, save.raw = T, file.name = 'test', parallel.sz=15) 
res_wide = make_plot_data(res)
p = plot_xCell(res_wide)

fig_fn = sprintf('%s/xCell_percentages.pdf',fig_dir)
save_plot(fig_fn,p,base_height=8,base_width=8)
