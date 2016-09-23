args=commandArgs(T)
eqtl_file=args[1]
gwas_file=args[2]
gene_id=args[3]
eqtl_file='../processed_data/160824/eCAVIAR_input/ENSG00000182511.7.eqtl.zscore'
gwas_file='../processed_data/160824/eCAVIAR_input/ENSG00000182511.7.gwas.zscore'
# gene_id='ENSG00000182511.7'

eqtl=fread(eqtl_file)%>%dplyr::rename(id=V1,zscore=V2)%>%dplyr::mutate(rank=rank(-abs(zscore),ties.method='random'))%>%as.data.table()
gwas=fread(gwas_file)%>%dplyr::rename(id=V1,zscore=V2)%>%dplyr::mutate(rank=rank(-abs(zscore),ties.method='random'))%>%as.data.table()

eqtl[rank<=100,]
gwas[rank<=100,]