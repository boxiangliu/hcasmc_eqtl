library(data.table)
library(dplyr)
library(stringr)
library(cowplot)

in_dir='../processed_data/160805/metasoft_output_subsample_52/'
fig_dir='../figures/160805/quality_control/'
if (!dir.exists(fig_dir)) dir.create(fig_dir,recursive=TRUE)
metasoft=fread(sprintf('%s/%s',in_dir,'metasoft_output.22.mcmc.txt'))
in_file=sprintf('%s/%s',in_dir,'metasoft_output.22.mcmc.txt')

read_metasoft=function(in_file){

	# read metasoft result: 
	header=scan(in_file,what='character',nlines=1,sep='\t')
	metasoft=fread(in_file,skip=1)


	# remove the extra column (issue due to white space)
	metasoft[,V107:=NULL]


	# read tissue name (study name):
	study_name=unlist(fread('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/Metasoft_tissue_order.alphabetical.txt',header=F))


	# set metasoft column names:
	col_names=c(paste('pvalue',study_name,sep='_'),paste('mvalue',study_name,sep='_'))
	col_names=c(header[1:16],col_names)
	stopifnot(ncol(metasoft)==length(col_names))
	setnames(metasoft,col_names)

	return(metasoft)
}

get_mvalue=function(metasoft){
	# subset to mvalues: 
	mvalue=metasoft%>%dplyr::select(contains('mvalue'))%>%as.data.frame()
	rownames(mvalue)=metasoft$RSID
	colnames(mvalue)=str_replace(colnames(mvalue),'mvalue_','')
	return(mvalue)
}

get_pvalue=function(metasoft){
	# subset to pvalues: 
	pvalue=metasoft%>%dplyr::select(contains('pvalue'))%>%dplyr::select(-c(1:5))%>%as.data.frame()
	rownames(pvalue)=metasoft$RSID
	colnames(pvalue)=str_replace(colnames(mvalue),'pvalue_','')
	return(pvalue)
}

get_tissue_specific_eqtl=function(mvalue,tissue){
	# select rows where HCASMC has m-value >= 0.9:
	is_eqtl=mvalue[,tissue]>=0.9


	# select rows where GTEx tissues all have m-values < 0.9:
	col_idx=which(colnames(mvalue)==tissue)
	max_mval_besides_this_tissue=apply(mvalue[,-col_idx],1,max,na.rm=T)
	is_not_eqtl_besides_this_tissue=max_mval_besides_this_tissue<0.9


	# select hcasmc-specific eQTLs:
	tissue_specific_eqtl_idx=which(is_eqtl&is_not_eqtl_besides_this_tissue)

	return(tissue_specific_eqtl_idx)
}

metasoft=read_metasoft(in_file)
mvalue=get_mvalue(metasoft)
pvalue=get_pvalue(metasoft)


# Plot number of tested genes per tissue:
container=list()
for (col in colnames(mvalue)){
	n_tested=sum(!is.na(mvalue[,col]))
	container[[col]]=data.frame(sample=col,n_tested)
}
n_tested=Reduce(rbind,container)
setDT(n_tested)
setorder(n_tested,n_tested)
n_tested[,sample:=factor(sample,levels=sample)]
p1=ggplot(n_tested,aes(sample,n_tested))+geom_bar(stat='identity')+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+ggtitle('Number of tested pairs')


# Subset to genes tested in HCASMC: 
mvalue2=mvalue[!is.na(mvalue$HCASMC),]
container=list()
for (col in colnames(mvalue2)){
	n_tested=sum(!is.na(mvalue2[,col]))
	container[[col]]=data.frame(sample=col,n_tested)
}
n_tested=Reduce(rbind,container)
setDT(n_tested)
setorder(n_tested,n_tested)
n_tested[,sample:=factor(sample,levels=sample)]
p2=ggplot(n_tested,aes(sample,n_tested))+geom_bar(stat='identity')+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+ggtitle('Number of tested pairs (subsetted to tested pairs in HCASMC)')


# Subset to genes tested in at least one other tissue besides HCASMC:
container=list()
for (tissue in colnames(mvalue2)){
	print(tissue)
	n_specific=length(get_tissue_specific_eqtl(mvalue2,tissue))
	container[[tissue]]=data.frame(tissue,n_specific)
}
n_specific=Reduce(rbind,container)
setDT(n_specific)
setorder(n_specific,n_specific)
n_specific[,tissue:=factor(tissue,levels=tissue)]
p3=ggplot(n_specific,aes(tissue,n_specific))+geom_bar(stat='identity')+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+ggtitle('Number of tissue specific pairs')


# Save plots:
pdf(sprintf('%s/quality_control.pdf',fig_dir))
p1;p2;p3
dev.off()