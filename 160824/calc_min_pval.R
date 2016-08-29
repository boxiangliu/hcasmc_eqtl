args=commandArgs(T)
eqtl_file=args[1]
gwas_file=args[2]
gene_id=args[3]
# eqtl_file='../processed_data/160824/eCAVIAR_input/ENSG00000182511.7.eqtl.zscore'
# gwas_file='../processed_data/160824/eCAVIAR_input/ENSG00000182511.7.gwas.zscore'
# gene_id='ENSG00000182511.7'

eqtl=fread(eqtl_file)%>%rename(id=V1,zscore=V2)%>%mutate(pval=2*(1-pnorm(abs(zscore))))
gwas=fread(gwas_file)%>%rename(id=V1,zscore=V2)%>%mutate(pval=2*(1-pnorm(abs(zscore))))
eqtl_min_pval=eqtl%>%summarize(min_pval=min(pval))
gwas_min_pval=gwas%>%summarize(min_pval=min(pval))


writeLines(paste(gene_id,gwas_min_pval$min_pval[1],eqtl_min_pval$min_pval[1],sep='\t'))
