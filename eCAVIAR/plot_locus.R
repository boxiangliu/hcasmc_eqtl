args=commandArgs(T)
gwas_file=args[1]
eqtl_file=args[2]
out_figure=args[3]
# gwas_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/eCAVIAR/eCAVIAR_input/ENSG00000100014.15.gwas.zscore'
# eqtl_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data//eCAVIAR/eCAVIAR_input/ENSG00000100014.15.eqtl.zscore'
# out_figure='../figures/eCAVIAR/ENSG00000100014.15.pdf'


# read gwas and eqtl data: 
gwas=fread(gwas_file)
setnames(gwas,c('id','zscore'))
eqtl=fread(eqtl_file)
setnames(eqtl,c('id','zscore'))

# calculate p-value: 
gwas[,pval_gwas:=-log10(2*(1-pnorm(abs(zscore))))]
eqtl[,pval_eqtl:=-log10(2*(1-pnorm(abs(zscore))))]

# merge gwas and eqtl: 
merge=merge(gwas,eqtl,by='id')%>%select(id,pval_gwas,pval_eqtl)

# melt: 
melt=melt(merge,id.vars='id',measure.vars=c('pval_gwas','pval_eqtl'),variable.name='type',value.name='pval')

# plot locus:
p=ggplot(melt,aes(x=id,y=pval,color=type))+geom_point()+theme_bw()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
save_plot(out_figure,p,base_height=6)