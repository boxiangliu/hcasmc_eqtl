library(stringr)
library(data.table)
library(dplyr)
library(dtplyr)
library(cowplot)

# Variables:
out_dir='../processed_data/gwas_eqtl_overlap/overlap.metasoft/'
fig_dir='../figures/gwas_eqtl_overlap/gwas_threshold/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

# read GWAS: 
gwas=fread('../data/gwas/CARDIoGRAMplusC4D/cad.add.160614.website.txt')
gwas=gwas%>%select(chrom=chr,pos=bp_hg19,rsid=markername,pval=p_dgc)


# constants: 
in_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/metasoft_output_subsample_52/metasoft_output.1.mcmc.txt'
N_TISSUE=45


# read metasoft result: 
header=scan('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/metasoft_output_subsample_52/metasoft_output.1.mcmc.txt',what='character',nlines=1,sep='\t')
metasoft=list()
for (i in seq(1,22)){
	in_file=sprintf('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/metasoft_output_subsample_52/metasoft_output.%s.mcmc.txt',i)
	print(sprintf('INFO - %s', in_file))
	metasoft[[i]]=fread(in_file,skip=1)
}
metasoft2=Reduce(rbind,metasoft)
metasoft=metasoft2

# remove the extra column (issue due to white space)
metasoft[,V107:=NULL]


# read tissue name (study name):
study_name=unlist(fread('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/Metasoft_tissue_order.alphabetical.txt',header=F))


# set metasoft column names:
col_names=c(paste('pvalue',study_name,sep='_'),paste('mvalue',study_name,sep='_'))
col_names=c(header[1:16],col_names)
stopifnot(ncol(metasoft)==length(col_names))
setnames(metasoft,col_names)
saveRDS(metasoft,sprintf('%s/metasoft.rds',out_dir))

# subset to mvalues: 
mvalue=metasoft%>%dplyr::select(contains('mvalue'))%>%as.data.frame()
stopifnot(ncol(mvalue)==N_TISSUE)
rownames(mvalue)=metasoft$RSID
colnames(mvalue)=str_replace(colnames(mvalue),'mvalue_','')


# subset to association with at least 10 tissues: 
mvalue_bak=mvalue
mvalue=mvalue_bak[metasoft$`#STUDY`>=10,]


# append pid, chrom, and pos: 
tmp=str_split_fixed(row.names(mvalue),'_',n=6)[,1:3]
tmp=data.frame(tmp,stringsAsFactors=F)
tmp[,2]=as.integer(tmp[,2])
tmp[,3]=as.integer(tmp[,3])
colnames(tmp)=c('pid','chrom','pos')
mvalue=cbind(tmp,mvalue)


# merge gwas and eQTL:
mvalue=as.data.table(mvalue)
merge=merge(gwas,mvalue,by=c('chrom','pos'))


# calculate overlap summary:
mvalue_threshold=0.9
overlap_summary=data.frame()
for (tissue in study_name){
	eqtl=merge[,tissue,with=F]
	sig_eqtl=eqtl>=mvalue_threshold
	n_sig_eqtl=sum(sig_eqtl,na.rm=T)
	for (gwas_threshold in c(1e-3,1e-4,1e-5,1e-6,1e-7)){
		sig_gwas=merge$pval<gwas_threshold
		n_sig_overlap=sum(sig_eqtl*sig_gwas,na.rm=T)
		tmp=data.frame(gwas_threshold=gwas_threshold,n_sig_eqtl=n_sig_eqtl,n_sig_overlap=n_sig_overlap,tissue=tissue)
		overlap_summary=rbind(overlap_summary,tmp)
	}
}


# calculate fraction of overlap over eQTLs:
overlap_summary=as.data.table(overlap_summary)
overlap_summary[,fraction:=n_sig_overlap/n_sig_eqtl]
overlap_summary[,tissue:=as.character(tissue)]


# Make plot:
pdf(sprintf('%s/overlap.pdf',fig_dir))
# fraction of overlap over eQTL:
ggplot(overlap_summary,aes(x=gwas_threshold,y=fraction,color=tissue,size=ifelse(tissue=='HCASMC'|tissue=='Artery_Coronary',2,1)))+geom_point()+scale_x_log10()+theme(axis.text.x=element_text(angle=45,hjust=1))+scale_color_discrete(guide=F)+scale_size(guide=F)+xlab('GWAS threshold')+ylab('Fraction of overlap')+geom_text(aes(label=ifelse(gwas_threshold==1e-3&(tissue=='HCASMC'|tissue=='Artery_Coronary'),tissue,'')),hjust=1.1)
ggplot(overlap_summary[gwas_threshold==1e-6],aes(y=reorder(tissue,fraction,FUN=mean),x=fraction))+geom_point()+scale_color_discrete(guide=F)+ylab('Tissue')+xlab('Fraction of overlap')


# just number of overlaps:
ggplot(overlap_summary,aes(x=gwas_threshold,y=n_sig_overlap,color=tissue,size=ifelse(tissue=='HCASMC'|tissue=='Artery_Coronary',2,1)))+geom_point()+scale_x_log10()+theme(axis.text.x=element_text(angle=45,hjust=1))+scale_color_discrete(guide=F)+scale_size(guide=F)+xlab('GWAS threshold')+ylab('Fraction of overlap')+geom_text(aes(label=ifelse(gwas_threshold==1e-3&(tissue=='HCASMC'|tissue=='Artery_Coronary'),tissue,'')),hjust=1.1)
ggplot(overlap_summary[gwas_threshold==1e-6],aes(y=reorder(tissue,n_sig_overlap,FUN=mean),x=n_sig_overlap))+geom_point()+scale_color_discrete(guide=F)+ylab('Tissue')+xlab('Fraction of overlap')


# just number of eQTLs:
ggplot(overlap_summary,aes(x=gwas_threshold,y=n_sig_eqtl,color=tissue,size=ifelse(tissue=='HCASMC'|tissue=='Artery_Coronary',2,1)))+geom_point()+scale_x_log10()+theme(axis.text.x=element_text(angle=45,hjust=1))+scale_color_discrete(guide=F)+scale_size(guide=F)+xlab('GWAS threshold')+ylab('Fraction of overlap')+geom_text(aes(label=ifelse(gwas_threshold==1e-3&(tissue=='HCASMC'|tissue=='Artery_Coronary'),tissue,'')),hjust=1.1)
ggplot(overlap_summary[gwas_threshold==1e-6],aes(y=reorder(tissue,n_sig_eqtl,FUN=mean),x=n_sig_eqtl))+geom_point()+scale_color_discrete(guide=F)+ylab('Tissue')+xlab('Fraction of overlap')
dev.off()


#-------------- Subset to genes tested in HCASMC -----------# 
# Subset to genes tested in HCASMC: 
mvalue=mvalue[!is.na(mvalue$HCASMC),]

# append pid, chrom, and pos: 
tmp=str_split_fixed(rownames(mvalue),'_',n=6)[,1:3]
tmp=data.frame(tmp,stringsAsFactors=F)
tmp[,2]=as.integer(tmp[,2])
tmp[,3]=as.integer(tmp[,3])
colnames(tmp)=c('pid','chrom','pos')
mvalue=cbind(tmp,mvalue)


# merge gwas and eQTL:
setDT(mvalue)
merge=merge(gwas,mvalue,by=c('chrom','pos'))


# Count the number of overlaps: 
container=list()
mvalue_threshold=0.9
i=0
for (gwas_threshold in c(1e-3,1e-4,1e-5,1e-6,1e-7)){
	sig_gwas=merge[pval<gwas_threshold,]

	for (tissue in study_name){
		i=i+1
		setnames(sig_gwas,tissue,'mvalue')
		max_mvalue=suppressWarnings(sig_gwas[,list(mvalue=max(mvalue,na.rm=TRUE)),by=c('chrom','pos')]) # warning can be safely ignored.
		max_mvalue=max_mvalue[!is.infinite(mvalue)]
		n_sig=max_mvalue[,sum(mvalue>=mvalue_threshold,na.rm=TRUE)]
		p=n_sig/nrow(sig_gwas)
		container[[i]]=data.frame(n=n_sig,p=p,tissue=tissue,gwas_threshold=gwas_threshold)
		setnames(sig_gwas,'mvalue',tissue)
	}
}
overlap_summary=Reduce(rbind,container)


pdf(sprintf('%s/overlap.genes_tested_in_hcasmc.pdf',fig_dir))
ggplot(overlap_summary,aes(gwas_threshold,p,color=tissue,alpha=ifelse(tissue=='HCASMC',10,0.1)))+geom_line()+scale_color_discrete(guide='none')+scale_x_log10(breaks=10^-seq(3,7))+scale_alpha_continuous(guide='none')+xlab('GWAS threshold')+ylab('% GWAS explained by eQTL')
dev.off()


