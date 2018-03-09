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

ukbb_fn='../data/gwas/ukbb/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz'
ukbb_marker_col='snptestid'
ukbb_pval_col='p-value_gc'
ukbb=read_full_metal(
	in_fn = ukbb_fn,
	marker_col = ukbb_marker_col,
	pval_col = ukbb_pval_col,
	a1_col = 'effect_allele',
	a2_col = 'noneffect_allele',
	effect_col = 'logOR',
	se_col = 'se_gc')


tcf21_fn='../processed_data/rasqual/output_pval/chr6/ENSG00000118526.6_TCF21.pval.txt'
tcf21 = read_full_metal(
	in_fn = tcf21_fn,
	marker_col = 'rsid',
	pval_col = 'pval',
	a1_col = 'ref',
	a2_col = 'alt',
	effect_col = 'pi',
	se_col = 'chisq'
	)

merged = merge(ukbb, tcf21, by = 'rsid')

x = ukbb[,list(rsid,a1,a2)]
y = tcf21[,list(rsid,a1,a2)]
align_genotype = function(x,y){
	z = merge(x,y,by='rsid')
	zx = paste(z$a1.x,z$a2.x,sep=':')
	zy = paste(z$a1.y,z$a2.y,sep=':')
	zy_flip = paste(z$a2.y,z$a1.y,sep=':')
	same = as.integer(zx == zy)
	opposite = as.integer(zx == zy_flip)


}