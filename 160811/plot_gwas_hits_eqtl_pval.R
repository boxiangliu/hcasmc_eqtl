# plot the -log10(pvalue) for each tested gene:snp pair at GWAS top hits: 

# read tested pairs: 
tested_pairs=fread('../processed_data/160811/tested_genes_at_gwas_top_hits.txt')%>%dplyr::select(1:6)
setnames(tested_pairs,c('pheno','geno','dist','pval','beta','se'))

# cast geno column to factors:
tested_pairs$geno=factor(tested_pairs$geno, levels=unique(tested_pairs$geno))


# make facted scatter plot: 
p=ggplot(tested_pairs,aes(x=pheno,y=-log10(pval)))+geom_point()+facet_grid(.~geno,scales='free')+theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+geom_hline(yintercept=-log10(0.05),color='red',linetype='dashed')
save_plot('../figures/160811/gwas_hits_gene_pval.pdf',p,base_width=200,base_height=5,limitsize=FALSE)