# read gwas data:
gwas=fread('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/gwas/cad.add.160614.website.txt')
gwas=gwas%>%mutate(gwas_zscore=beta/se_dgc)%>%dplyr::select(chr,pos=bp_hg19,effect_allele,noneffect_allele,gwas_zscore,markername)
gwas=as.data.table(gwas)


# read eqtl file names:
in_files=list.files(path='../processed_data/160824/eqtl_by_gene',pattern='eqtl.zscore',full.names=T)
for (in_file in in_files){
	# get gene id:
	gene_id=basename(in_file)%>%str_replace('.eqtl.zscore','')
	message(gene_id	)

	# reformat eqtl: 
	eqtl=fread(in_file)
	setnames(eqtl,c('id','eqtl_zscore'))
	tmp=as.data.table(str_split_fixed(eqtl$id,"_",5)[,1:4])
	setnames(tmp,c('chr','pos','ref','alt'))
	tmp$chr=as.integer(tmp$chr)
	tmp$pos=as.integer(tmp$pos)
	eqtl=cbind(eqtl,tmp)


	# intersect gwas and eqtl:
	gwas_small=gwas[chr==unique(eqtl$chr)&pos>=min(eqtl$pos)&pos<=max(eqtl$pos),]
	merge=merge(eqtl,gwas_small,by=c('chr','pos'))
	merge[,ref2:=ifelse(nchar(ref)!=nchar(alt),'I',ref)]
	merge[,alt2:=ifelse(nchar(ref)!=nchar(alt),'D',alt)]
	merge=merge[which(pmin(merge$ref2,merge$alt2)==pmin(merge$effect_allele,merge$noneffect_allele)),]
	

	# write output:
	write.table(merge[,.(id,eqtl_zscore)],paste0('../processed_data/eCAVIAR_input/',gene_id,'.eqtl.zscore'),quote=F,row.names=F,col.names=F,sep='\t')
	write.table(merge[,.(id,gwas_zscore)],paste0('../processed_data/eCAVIAR_input/',gene_id,'.gwas.zscore'),quote=F,row.names=F,col.names=F,sep='\t')
	write.table(merge[,markername],paste0('../processed_data/eCAVIAR_input/',gene_id,'.markername'),quote=F,row.names=F,col.names=F,sep='\t')
}
