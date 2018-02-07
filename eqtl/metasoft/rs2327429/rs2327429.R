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
out_dir='../processed_data/eqtl/metasoft/rs2327429/rs2327429/'
if (!dir.exists(out_dir)){dir.create(out_dir,recursive=TRUE)}

# Functions:
read_eqtl=function(in_fn){
	eqtl=fread(sprintf('gunzip -c %s | grep ENSG00000118526.6 | grep 6_134209837_T_C_b37',in_fn),col.names=c('fid','sid','dist','pval','beta','se'))
	return(eqtl)
}

select_eSNP=function(eqtl,subset){
	eqtl_subset=eqtl[fid==subset$fid&sid==subset$sid,]
	return(eqtl_subset)
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
	fwrite(merged,paste0(out_dir,'/',tissue,'.rs2327429.txt'),sep='\t',col.names=FALSE)
}

main=function(){
	eqtl=read_eqtl(in_fn)
	top_eSNP=select_eSNP(eqtl,subset=list(fid='ENSG00000118526.6',sid='6_134209837_T_C_b37'))
	top_eSNP$tissue='HCASMC'
	fwrite(top_eSNP,paste0(out_dir,'/HCASMC.rs2327429.txt'),sep='\t',col.names=FALSE)

	gtex_tissue_dirs=get_gtex_tissue_dirs(gtex_dir)
	foreach(i=1:length(gtex_tissue_dirs))%dopar%{
		message(gtex_tissue_dirs[i])
		subset_gtex_eqtls(gtex_tissue_dirs[i],subset=top_eSNP[,list(fid,sid)],out_dir)
	}
}

main()