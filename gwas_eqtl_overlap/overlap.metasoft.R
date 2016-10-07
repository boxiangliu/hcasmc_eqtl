library(stringr)
library(data.table)
library(dplyr)
library(dtplyr)
library(cowplot)

# read GWAS: 
gwas=fread('../data/gwas/cad.add.160614.website.txt')


# format GWAS data:
gwas=gwas%>%select(chrom=chr,pos=bp_hg19,rsid=markername,pval=p_dgc)



in_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/metasoft_output_subsample_52/metasoft_output.1.mcmc.txt'
# constants: 
N_TISSUE=45


# read metasoft result: 
header=scan('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/metasoft_output_subsample_52/metasoft_output.1.mcmc.txt',what='character',nlines=1,sep='\t')
metasoft=data.table()
for (i in seq(1,22)){
	in_file=sprintf('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/metasoft_output_subsample_52/metasoft_output.%s.mcmc.txt',i)
	metasoft=rbind(metasoft,fread(in_file,skip=1))
}


# remove the extra column (issue due to white space)
metasoft[,V107:=NULL]


# read tissue name (study name):
study_name=unlist(fread('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/Metasoft_tissue_order.alphabetical.txt',header=F))


# set metasoft column names:
col_names=c(paste('pvalue',study_name,sep='_'),paste('mvalue',study_name,sep='_'))
col_names=c(header[1:16],col_names)
stopifnot(ncol(metasoft)==length(col_names))
setnames(metasoft,col_names)


# subset to pvalues:
pvalue=metasoft%>%dplyr::select(contains('pvalue'))%>%dplyr::select(-c(1:5))%>%as.data.frame()
rownames(pvalue)=metasoft$RSID
colnames(pvalue)=str_replace(colnames(pvalue),'pvalue_','')


# subset to association with at least 10 tissues: 
pvalue_bak=pvalue
pvalue=pvalue_bak[metasoft$`#STUDY`>=10,]


# append pid and sid: 
tmp=str_split_fixed(row.names(pvalue),'_',n=6)[,1:3]
tmp=data.frame(tmp,stringsAsFactors=F)
tmp[,2]=as.integer(tmp[,2])
tmp[,3]=as.integer(tmp[,3])
colnames(tmp)=c('pid','chrom','pos')
pvalue_bak=pvalue
pvalue=cbind(tmp,pvalue_bak)

tmp

# merge gwas and eQTL:
pvalue=as.data.table(pvalue)
merge=merge(gwas,pvalue,by=c('chrom','pos'))


overlap_summary=data.frame()
for (eqtl_threshold in c(1e-3,1e-4,1e-5,1e-6,1e-7)){
	for (tissue in study_name){
		eqtl=merge[,tissue,with=F]
		sig_eqtl=eqtl<eqtl_threshold
		n_sig_eqtl=sum(sig_eqtl,na.rm=T)
		for (gwas_threshold in c(1e-3,1e-4,1e-5,1e-6,1e-7)){
			sig_gwas=merge$pval<gwas_threshold
			n_sig_overlap=sum(sig_eqtl*sig_gwas,na.rm=T)
			tmp=data.frame(gwas_threshold=gwas_threshold,eqtl_threshold=eqtl_threshold,n_sig_eqtl=n_sig_eqtl,n_sig_overlap=n_sig_overlap,tissue=tissue)
			overlap_summary=rbind(overlap_summary,tmp)
		}
	}
}


overlap_summary=as.data.table(overlap_summary)
overlap_summary[,fraction:=n_sig_overlap/n_sig_eqtl]
ggplot(overlap_summary[eqtl_threshold==1e-5,],aes(x=gwas_threshold,y=fraction,color=tissue,size=ifelse(tissue=='HCASMC',2,1)))+geom_point()+scale_x_log10()+theme(axis.text.x=element_text(angle=45,hjust=1))+scale_color_discrete(guide=F)



# subset to mvalues: 
mvalue=metasoft%>%dplyr::select(contains('mvalue'))%>%as.data.frame()
stopifnot(ncol(mvalue)==N_TISSUE)
rownames(mvalue)=metasoft$RSID
colnames(mvalue)=str_replace(colnames(mvalue),'mvalue_','')


# subset to association with at least 10 tissues: 
mvalue_bak=mvalue
mvalue=mvalue_bak[metasoft$`#STUDY`>=10,]


# append pid and sid: 
tmp=str_split_fixed(row.names(mvalue),'_',n=6)[,1:3]
tmp=data.frame(tmp,stringsAsFactors=F)
tmp[,2]=as.integer(tmp[,2])
tmp[,3]=as.integer(tmp[,3])
colnames(tmp)=c('pid','chrom','pos')
mvalue_bak=mvalue
mvalue=cbind(tmp,mvalue_bak)


# merge gwas and eQTL:
mvalue=as.data.table(mvalue)
merge=merge(gwas,mvalue,by=c('chrom','pos'))

overlap_summary=data.frame()
for (tissue in study_name){
	eqtl=merge[,tissue,with=F]
	sig_eqtl=eqtl>=0.9
	n_sig_eqtl=sum(sig_eqtl,na.rm=T)
	for (gwas_threshold in c(1e-3,1e-4,1e-5,1e-6,1e-7)){
		sig_gwas=merge$pval<gwas_threshold
		n_sig_overlap=sum(sig_eqtl*sig_gwas,na.rm=T)
		tmp=data.frame(gwas_threshold=gwas_threshold,n_sig_eqtl=n_sig_eqtl,n_sig_overlap=n_sig_overlap,tissue=tissue)
		overlap_summary=rbind(overlap_summary,tmp)
	}
}


overlap_summary=as.data.table(overlap_summary)
overlap_summary[,fraction:=n_sig_overlap/n_sig_eqtl]

ggplot(overlap_summary,aes(x=gwas_threshold,y=fraction,color=tissue,size=ifelse(tissue=='HCASMC',2,1)))+geom_point()+scale_x_log10()+theme(axis.text.x=element_text(angle=45,hjust=1))+scale_color_discrete(guide=F)
ggplot(overlap_summary[gwas_threshold==1e-3],aes(x=reorder(tissue,fraction,FUN=mean),y=fraction))+geom_point()+theme(axis.text.x=element_text(angle=45,hjust=1))+scale_color_discrete(guide=F)
save.image('../processed_data/gwas_eqtl_overlap/overlap.metasoft.RData')
