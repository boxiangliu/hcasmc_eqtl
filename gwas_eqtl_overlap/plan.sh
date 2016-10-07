# read gwas data
# read eqtl data hcasmc and some related tissues
# for gwas pvalue threshold: 
# 	filter gwas
# 	for each eqtl:
# 		overlap gwas with eqtl
#		take the most significant eqtl
# 		make qqplot 
# 		append to overlap summary


# overlap.metasoft.R 
# read gwas data
# for each chromosome: 
# 	read metasoft for that chromosome
# 	extract pvalue 
# 	subset to associations with over 10 tissues(?) 
#	extract the sid and pid
# 	extract chr and pos
#	overlap the eQTL and GWAS
# 	for tissue
# 		extract tissue -> eqtl
#		sig_eqtl=eqtl<eqtl_threshold
#		n_sig_eqtl=sum(sig_eqtl)
#		for gwas_threshold: 
#			sig_gwas=gwas_pval<gwas_threshold
# 			n_sig_overlap=sum(sig_eqtl*sig_gwas)
# 			append to overlap summary




