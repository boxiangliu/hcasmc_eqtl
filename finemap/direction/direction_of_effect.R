# Determine direction of effect
# Boxiang Liu
# 2018-03-08
library(data.table)
# library(devtools)
# detach('package:locuscomparer',unload=TRUE)
# devtools::install_github('boxiangliu/locuscomparer')
library(locuscomparer)
library(cowplot)
library(ggrepel)
library(stringr)

# Variables:
PVAL = 0.001
ukbb_fn='../data/gwas/ukbb/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz'
ukbb_marker_col='snptestid'
ukbb_pval_col='p-value_gc'
fig_dir = '../figures/finemap/direction/direction_of_effect/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

eqtl_fn = c(
	TCF21 = '../processed_data/rasqual/output_pval/chr6/ENSG00000118526.6_TCF21.pval.txt',
	SMAD3 = '../processed_data/rasqual/output_pval/chr15/ENSG00000166949.11_SMAD3.pval.txt',
	FES = '../processed_data/rasqual/output_pval/chr15/ENSG00000182511.7_FES.pval.txt',
	PDGFRA = '../processed_data/rasqual/output_pval/chr4/ENSG00000134853.7_PDGFRA.pval.txt',
	SIPA1 = '../processed_data/rasqual/output_pval/chr11/ENSG00000213445.4_SIPA1.pval.txt')

# Functions:
merge_and_align_genotype = function(x,y,y_flip_effect = function(b) {1 - b}){
	shared_rsid = merge(x[,list(rsid)],y[,list(rsid)],by='rsid')$rsid
	x_sub = x[rsid %in% shared_rsid]
	y_sub = y[rsid %in% shared_rsid]
	same = merge(x_sub,y_sub,by=c('rsid','a1','a2'))
	message = sprintf('%s SNPs are in the same direction',nrow(same))
	message(message)

	setnames(y_sub,c('a1','a2'),c('a2','a1'))
	y_sub[,effect:=y_flip_effect(effect)]
	flipped = merge(x_sub,y_sub,by=c('rsid','a1','a2'))
	message = sprintf('%s SNPs are in the opposite direction',nrow(flipped))
	message(message)

	merged = rbind(same,flipped)
	message = sprintf('%s SNPs are dropped',length(shared_rsid) - nrow(merged))
	message(message)

	return(merged)
}

plot_effect_sizes = function(merged,gene_name){
	ylabel = sprintf('%s\nAllelic ratio',gene_name)
	ggplot(merged[pval.x<PVAL&pval.y<PVAL],aes(effect.x,effect.y))+
		geom_point(size = 3, alpha = 0.5) + 
		stat_smooth(method='lm', formula = I(y-0.5) ~ x + 0, position = position_nudge(y = 0.5)) + 
		xlab('UKBB log(OR)') + 
		ylab(ylabel)
}

plot_blank = function(){
	ggplot()+geom_blank(aes(1,1))+
		theme(
		plot.background = element_blank(), 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		panel.border = element_blank(),
		panel.background = element_blank(),
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		axis.text.x = element_blank(), 
		axis.text.y = element_blank(),
		axis.ticks = element_blank(),
		axis.line = element_blank()
	)
}
# Main: 
ukbb=read_full_metal(
	in_fn = ukbb_fn,
	marker_col = ukbb_marker_col,
	pval_col = ukbb_pval_col,
	a1_col = 'effect_allele',
	a2_col = 'noneffect_allele',
	effect_col = 'logOR',
	se_col = 'se_gc')


plot_list = list()
for (i in seq_along(eqtl_fn)){
	fn = eqtl_fn[i]
	gene_name = names(eqtl_fn)[i]
	message(gene_name)

	eqtl = read_full_metal(
		in_fn = fn,
		marker_col = 'rsid',
		pval_col = 'pval',
		a1_col = 'ref',
		a2_col = 'alt',
		effect_col = 'pi',
		se_col = 'chisq'
		)
	
	merged = merge_and_align_genotype(ukbb,eqtl,y_flip_effect=function(b) {1 - b})

	p = plot_effect_sizes(merged,gene_name)
	plot_list[[i]] = p
	fig_fn = sprintf('%s/%s_effect_size.pdf',fig_dir,gene_name)
	save_plot(fig_fn,p)
}


row1 = plot_grid(plot_list[[1]],plot_list[[2]],labels = c('A','B'),align = 'h')
row2 = plot_grid(plot_list[[3]],plot_list[[4]],labels = c('C','D'),align = 'h')
row3 = plot_grid(plot_list[[5]],plot_blank(),labels = c('E',''))

multi_panel = plot_grid(row1,row2,row3,labels='',ncol=1)
fig_fn = sprintf('%s/multi_panel_effect_size.pdf',fig_dir)
save_plot(fig_fn,multi_panel,base_width=7.5,base_height=9)