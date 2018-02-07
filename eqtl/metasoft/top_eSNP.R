# Select top eSNP per gene
# Boxiang Liu
# 2018-01-05

library(data.table)
library(foreach)
library(doMC)
registerDoMC(5)

# Variables:
in_fn='../processed_data/eqtl/fastqtl/output/nominal/all.txt.gz'
gtex_dir='/srv/persistent/bliu2/HCASMC_eQTL/processed_data//160816/subsampling/'
out_dir='../processed_data/eqtl/metasoft/top_eSNP/'
if (!dir.exists(out_dir)){dir.create(out_dir,recursive=TRUE)}

# Functions:
read_eqtl=function(in_fn){
	eqtl=fread(sprintf('gunzip -c %s',in_fn),col.names=c('fid','sid','dist','pval','beta','se'))
	return(eqtl)
}

select_top_eSNP=function(eqtl){
	eqtl[,min_pval:=min(pval),by='fid']
	top_eqtl=eqtl[pval==min_pval,list(fid,sid,dist,pval,beta,se)]
	return(top_eqtl)
}

get_gtex_tissue_dirs=function(gtex_dir){
	tmp=list.dirs(gtex_dir)
	gtex_tissue_dirs=tmp[tmp!=gtex_dir]
	return(gtex_tissue_dirs)
}

subset_gtex_eqtls=function(in_dir,subset,out_dir,tissue=basename(in_dir)){
	print(paste0('INFO - tissue:', tissue))

	in_fn=list.files(in_dir,pattern='allpairs.txt.gz',full.names=TRUE)
	stopifnot(length(in_fn)==1)
	print(paste0('INFO - input: ',in_fn))

	eqtl=read_eqtl(in_fn)
	merged=merge(eqtl,subset,by=c('fid','sid'),all=FALSE)
	print(paste0('INFO - ',nrow(merged),' rows remaining.'))
	merged$tissue=tissue
	fwrite(merged,paste0(out_dir,'/',tissue,'.top_eSNP.txt'),sep='\t',col.names=FALSE)
}

main=function(){
	eqtl=read_eqtl(in_fn)
	top_eSNP=select_top_eSNP(eqtl)
	top_eSNP$tissue='HCASMC'
	fwrite(top_eSNP,paste0(out_dir,'/HCASMC.top_eSNP.txt'),sep='\t',col.names=FALSE)

	gtex_tissue_dirs=get_gtex_tissue_dirs(gtex_dir)
	foreach(i=1:length(gtex_tissue_dirs))%dopar%{
		subset_gtex_eqtls(gtex_tissue_dirs[i],subset=top_eSNP[,list(fid,sid)],out_dir)
	}
}

main()




