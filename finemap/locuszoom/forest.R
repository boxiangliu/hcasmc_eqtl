# Make forest plot
# Boxiang Liu
# 2017-12-07

library(stringr)
library(foreach)
library(data.table)
library(cowplot)

args=commandArgs(T)
gene_id=args[1]
snp_id=args[2]
in_dir=args[3]
fig_dir=args[4]

gene_id='ENSG00000166949.11'
snp_id='rs17293632'
in_dir='/srv/persistent/bliu2/HCASMC_eQTL/figures/finemap/locuszoom/'
fig_dir='/srv/persistent/bliu2/HCASMC_eQTL/figures/finemap/locuszoom/forest/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

get_gtex_eqtl=function(in_dir,gene_id,snp_id){
	eqtl_fn_list=list.files(in_dir,pattern='fastQTL.txt$',recursive=TRUE,full.names=TRUE)
	eqtl_fn_list=eqtl_fn_list[str_detect(eqtl_fn_list,gene_id)]

	eqtl=foreach(eqtl_fn=eqtl_fn_list,.combine='rbind')%do%{
		x=fread(eqtl_fn,col.names=c('gene_id','snp','dist','pval','beta','se'))
		y=x[snp==snp_id,list(gene_id,snp,pval,beta,se)]
		tissue=str_extract(eqtl_fn,'(?<=ENSG.{10,30}/)(.+?)(?=/fastQTL)')
		y$tissue=tissue
		return(y)
	}
	return(eqtl)
}


in_dir='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/rasqual/output_pval/'
get_hcasmc_eqtl=function(in_dir,gene_id,snp_id){
	eqtl_fn=list.files(in_dir,gene_id,recursive=TRUE,full.names=TRUE)
	eqtl=fread(eqtl_fn)
	eqtl=eqtl[rsid==snp_id]

	eqtl[,beta:=2*pi-1]
	eqtl[,z:=sign(beta)*abs(qnorm(pval/2))]
	eqtl[,se:=beta/z]

	return(eqtl[,list(gene_id=fid,snp=rsid,pval,beta,se,tissue='HCASMC')])
}


make_forest_plot=function(eqtl,fig_dir,gene_id,snp_id){
	setorder(eqtl,beta)
	eqtl[,tissue:=factor(tissue,levels=tissue)]
	pdf(sprintf('%s/%s_%s_forest.pdf',fig_dir,gene_id,snp_id))
	p=ggplot(eqtl,aes(y=beta,ymin=beta-ci,ymax=beta+ci,x=tissue))+geom_pointrange()+coord_flip()+geom_hline(yintercept=0,color='Red',linetype=2)
	print(p)
	dev.off()
}


make_pval_plot=function(eqtl,fig_dir,gene_id,snp_id){
	setorder(eqtl,pval)
	eqtl[,tissue:=factor(tissue,levels=tissue)]
	p=ggplot(eqtl,aes(x=-log10(pval),y=tissue))+geom_point()+geom_vline(xintercept=-log10(0.05/45),color='Red',linetype=2)+geom_vline(xintercept=-log10(0.05),color='pink',linetype=2)
	pdf(sprintf('%s/%s_%s_pval.pdf',fig_dir,gene_id,snp_id))
	print(p)
	dev.off()
}

gtex_eqtl=get_gtex_eqtl(in_dir,gene_id,snp_id)
hcasmc_eqtl=get_hcasmc_eqtl(in_dir,gene_id,snp_id)
eqtl=rbind(gtex_eqtl,hcasmc_eqtl)
eqtl[,ci:=1.96*se]



make_pval_plot(eqtl,fig_dir,gene_id,snp_id)