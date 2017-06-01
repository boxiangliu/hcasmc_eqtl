library(cowplot)
library(data.table)
library(dplyr)
library(dtplyr)
library(pryr)
library(gap)


# function to read eQTL files: 
readEQTL=function(command){
	eqtl=fread(command,select=c(1,2,3,7))
	message(nrow(eqtl),' associations.')
	setnames(eqtl,c('pid','chrom','pos','pval'))
	eqtl_minPval=eqtl%>%group_by(chrom,pos)%>%summarize(pval_min=min(pval))
	rm(eqtl)
	message(nrow(eqtl_minPval),' unique variants.')
	return(eqtl_minPval)
}

readEQTL2=function(command){
	eqtl=fread(command,select=c(1,2,3,7))
	message(nrow(eqtl),' associations.')
	setnames(eqtl,c('pid','chrom','pos','pval'))
	return(eqtl)
}


# read GWAS: 
gwas=fread('../data/gwas/cad.add.160614.website.txt')


# format GWAS data:
gwas=gwas%>%select(chrom=chr,pos=bp_hg19,rsid=markername,pval=p_dgc)


# read list of tissues: 
tissues=unlist(fread('../data/gtex/gtex.v6p.eqtl.tissues.with_hcasmc.txt',header=F))


# read eQTL data:
overlap_summary=data.frame()
# for (this_tissue in tissues){
for (this_tissue in c('HCASMC','Muscle_Skeletal')){
	
	command=sprintf('cat ../processed_data/160816/subsampling/%s/%s_52.allpairs.sid_parsed.*.txt',this_tissue,this_tissue)
	eqtl_minPval=readEQTL(command)
	for (thd1 in c(1e-3,1e-4,1e-5,1e-6,1e-7)){
		# thresholding gwas by pvalue:
		message('gwas threshold: ', thd1)
		gwas_thd=gwas[pval<thd1,]
		n_gwas=nrow(gwas_thd)

		# overlap eQTL and GWAS: 
		merged=merge(gwas_thd,eqtl_minPval,by=c('chrom','pos'))
		

		for (thd2 in c(1e-3,1e-4,1e-5,1e-6,1e-7)){
			message('eqtl threshold: ', thd2)
			# count significant eQTLs in the overlap: 
			n_overlap=sum(merged$pval_min<thd2)
			n_eqtl=sum(eqtl_minPval$pval_min<thd2)
			tmp_df=data.frame(n_eqtl=n_eqtl,n_overlap=n_overlap,n_gwas=n_gwas,gwas_thd=thd1,eqtl_thd=thd2,tissue=this_tissue)
			overlap_summary=rbind(overlap_summary,tmp_df)
		}
	}
	save(overlap_summary,file='../processed_data/gwas_eqtl_overlap/overlap.RData')
	rm(eqtl_minPval)
}

merged_hcasmc=merge(gwas,eqtl_hcasmc_minPval,by=c('chrom','pos'))
merged_muscle_skeletal=merge(gwas,eqtl_muscle_skeletal_minPval,by=c('chrom','pos'))
merged_muscle_skeletal
save(overlap_summary,file='../processed_data/gwas_eqtl_overlap/overlap.RData')
load('../processed_data/gwas_eqtl_overlap/overlap.RData')

overlap_summary=as.data.table(overlap_summary)
overlap_summary[,fraction:=n_overlap/n_eqtl]
ggplot(overlap_summary[eqtl_thd==1e-3&gwas_thd==1e-7,],aes(y=fraction,x=reorder(tissue,fraction,FUN=mean)))+geom_point()+theme(axis.text.x=element_text(angle=45,hjust=1))
overlap_summary[eqtl_thd==1e-3&gwas_thd==1e-7,c(1,2,6),with=F]


this_tissue='HCASMC'
command=sprintf('cat ../processed_data/160816/subsampling/%s/%s_52.allpairs.sid_parsed.*.txt',this_tissue,this_tissue)
eqtl_hcasmc=readEQTL2(command)
setnames(eqtl_hcasmc,'pval','pval_eqtl')
gwas_1e_6=gwas[pval<1e-6,]


merged_hcasmc_1e_6=merge(gwas_1e_6,eqtl_hcasmc,by=c('chrom','pos'))
qqunif(merged_hcasmc_1e_6$pval_eqtl)


this_tissue='Muscle_Skeletal'
command=sprintf('cat ../processed_data/160816/subsampling/%s/%s_52.allpairs.sid_parsed.*.txt',this_tissue,this_tissue)

eqtl_muscle_skeletal=readEQTL2(command)
setnames(eqtl_muscle_skeletal,'pval','pval_eqtl')

merged_muscle_skeletal_1e_6=merge(gwas_1e_6,eqtl_muscle_skeletal,by=c('chrom','pos'))
qq_muscle_skeletal=qqunif(merged_muscle_skeletal_1e_6$pval_eqtl,plot.it=F)
points(qq_muscle_skeletal,col='red')



eqtl_hcasmc_1e_6=eqtl_hcasmc[pval_eqtl<1e-6,]
eqtl_muscle_skeletal_1e_6=eqtl_muscle_skeletal[pval_eqtl<1e-6]
merged_hcasmc_1e_6_2=merge(gwas,eqtl_hcasmc_1e_6,by=c('chrom','pos'))
merged_muscle_skeletal_1e_6_2=merge(gwas,eqtl_muscle_skeletal_1e_6,by=c('chrom','pos'))

qqunif(merged_hcasmc_1e_6_2$pval,xlim=c(0,10),ylim=c(0,10))
qq_muscle_skeletal_2=qqunif(merged_muscle_skeletal_1e_6_2$pval,plot.it=F)
points(qq_muscle_skeletal,col='red')
ggplot(merged_hcasmc_1e_6_2,aes(pvals))+geom_density()+geom_density(data=merged_muscle_skeletal_1e_6,aes(pval))
hist(merged_muscle_skeletal_1e_6_2$pval)
hist(merged_hcasmc_1e_6_2$pval)
