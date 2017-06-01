library(data.table)
library(dplyr)
library(stringr)
library(doMC)
library(foreach)
library(cowplot)
registerDoMC(40)

atac_fn='../processed_data/hcasmc_specific_open_chromatin/peak_specificity_filt/HCASMC.bed'
metasoft_dir='../processed_data/160805/metasoft_output_subsample_52/'
fig_dir='../figures/eqtl_and_atacseq/overlap/specific_atacseq_and_metasoft/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}


# Function: 
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

get_tissue_specific_eqtl=function(mvalue,tissue,return_idx=TRUE){
	# select rows where HCASMC has m-value >= 0.9:
	is_eqtl=mvalue[,tissue]>=0.9


	# select rows where GTEx tissues all have m-values < 0.9:
	col_idx=which(colnames(mvalue)==tissue)
	max_mval_besides_this_tissue=apply(mvalue[,-col_idx],1,max,na.rm=T)
	is_not_eqtl_besides_this_tissue=max_mval_besides_this_tissue<0.9


	# select hcasmc-specific eQTLs:
	tissue_specific_eqtl_idx=which(is_eqtl&is_not_eqtl_besides_this_tissue)

	if (return_idx){
		return(tissue_specific_eqtl_idx)
	} else {
		specific_eqtl=data.table(id=names(tissue_specific_eqtl_idx),mvalue=mvalue2[tissue_specific_eqtl_idx,tissue])
		tmp=data.table(specific_eqtl[,str_split_fixed(id,'_',2)])
		tmp2=data.table(tmp[,str_split_fixed(V2,'_',5)[,1:4]])
		specific_eqtl=cbind(specific_eqtl,tmp,tmp2)
		setnames(specific_eqtl,c('id','mvalue','gene_id','sid','chr','pos','ref','alt'))
		return(specific_eqtl)
	}
}


# Read mvalue:  
container=list()
for (i in seq(22)){
	in_file=sprintf('%s/metasoft_output.%s.mcmc.txt',metasoft_dir,i)
	metasoft=read_metasoft(in_file)
	mvalue=get_mvalue(metasoft)
	mvalue2=mvalue[!is.na(mvalue$HCASMC),]
	container[[i]]=mvalue2
}
mvalue=Reduce(rbind,container)


# Get HCASMC mvalue rank:
set.seed(42)
cores=40
size=ceiling(nrow(mvalue)/cores)
mvalue_rank=foreach (i = 1:40,.combine='c') %dopar% {
	start=(i-1)*size+1
	end=min(i*size,nrow(mvalue))
	apply(mvalue[start:end,],1,function(x) {y=rank(-x, ties.method='random');y[25]})
}
id=data.table(id=names(mvalue_rank))
tmp=data.table(id[,str_split_fixed(id,'_',2)])
tmp2=data.table(tmp[,str_split_fixed(V2,'_',5)[,1:4]])
rank=cbind(id,tmp,tmp2,data.table(mvalue_rank))
setnames(rank,c('id','gene_id','sid','chr','pos','ref','alt','rank'))
rank[,chr:=paste0('chr',chr)]
rank[,c('start','end'):=as.integer(pos)]
setkey(rank,chr,start,end)


atac=fread(atac_fn)
atac[,c('start','end'):=list(start-500,end+500)]
setkey(atac,chr,start,end)
overlap=foverlaps(rank,atac,nomatch=0)


# Calculate mean rank:
mean_rank=overlap[,list(rank=mean(rank)),by='num_tissue']
setorder(mean_rank,num_tissue)
fit=lm(rank~num_tissue,data=mean_rank)
p=ggplot(mean_rank,aes(num_tissue,rank))+geom_point()+stat_smooth(method='lm')+annotate(geom='text',x=30,y=25,label=sprintf('y=%.03f+%.03ex\nP=%.03e',coefficients(fit)[1],coefficients(fit)[2],summary(fit)$coefficients[2,4]))+xlab('Number of Tissues Sharing the Peak')+ylab('HCASMC mvalue rank\n(lower=more specific)')


# Remove an outlier point and replot:
mean_rank2=mean_rank[-which(rank==max(rank)),]
fit=lm(rank~num_tissue,data=mean_rank2)
p2=ggplot(mean_rank2,aes(num_tissue,rank))+geom_point()+stat_smooth(method='lm')+annotate(geom='text',x=30,y=24.7,label=sprintf('y=%.03f+%.03ex\nP=%.03e',coefficients(fit)[1],coefficients(fit)[2],summary(fit)$coefficients[2,4]))+xlab('Number of Tissues Sharing the Peak')+ylab('HCASMC mvalue rank\n(lower=more specific)')


# Save plots:
pdf(sprintf('%s/rank_vs_num_tissues.pdf',fig_dir),width=4.4,height=4)
p;p2
dev.off()

