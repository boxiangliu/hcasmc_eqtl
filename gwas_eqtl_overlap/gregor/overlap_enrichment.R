library(data.table)
library(sinib)
library(dplyr)
library(stringr)
library(cowplot)
library(gplots)

# Variables: 
metasoft_dir='../processed_data/160805/metasoft_output_subsample_52/'
fig_dir='../figures/gwas_eqtl_overlap/gregor/overlap_enrichment/'
if (!dir.exists(fig_dir)) dir.create(fig_dir)

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


# Read GWAS and matched background SNPs (and LD SNPs):
ld_set=fread('../processed_data/gwas_atacseq_overlap/tmp/ld_set.tsv')
ld_set[,c('start','end'):=list(pos,pos)]
setkey(ld_set,chr,start,end)


# Read metasoft:
container=list()
for (i in seq(22)){
	in_file=sprintf('%s/metasoft_output.%s.mcmc.txt',metasoft_dir,i)
	metasoft=read_metasoft(in_file)
	mvalue=get_mvalue(metasoft)
	container[[i]]=mvalue
}
mvalue=Reduce(rbind,container)


# Caculate overlap enrichment:
id=data.table(id=rownames(mvalue))
rowdata=data.table(id[,str_split_fixed(id,'_',6)][,1:3])
setnames(rowdata,c('gene_id','chr','pos'))
rowdata[,chr:=paste0('chr',chr)]
rowdata[,c('start','end'):=as.integer(pos)]
rowdata$id=id

setkey(ld_set,chr,start,end)
container=list()
for (i in 1:ncol(mvalue)){
	print(sprintf('INFO - %s',colnames(mvalue)[i]))
	x=rowdata[which(mvalue[,i]>=0.9),]
	y=unique(x[,list(chr,start,end)])
	setkey(y,chr,start,end)

	# Overlap: 
	overlap=unique(foverlaps(ld_set,y))
	overlap[,c('i.start','i.end'):=NULL]


	overlap[,snp_overlap:=!is.na(start)]
	overlap[,loci_overlap:=any(snp_overlap),by='loci_index']
	overlap[,c('start','end'):=NULL]
	overlap=unique(overlap)
	stopifnot(nrow(overlap)==nrow(ld_set))


	overlap=overlap[ld_proxy==FALSE,]
	overlap[,p:=mean(loci_overlap),by='gwas_index']

	p=overlap[,list(p=unique(p)),by='gwas_index']
	n=rep(1,length(p$p))
	s=sum(overlap[snpID==gwas_index,loci_overlap])
	container[[i]]=data.table(tissue=colnames(mvalue)[i],pval=psinib(q=as.integer(s-1),size=as.integer(n),prob=p$p,lower.tail=F))
}
pval=Reduce(rbind,container)
setorder(pval,pval)
pval[,tissue:=factor(tissue,levels=tissue)]
p1=ggplot(pval,aes(tissue,-log10(pval)))+geom_point()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
save_plot(sprintf('%s/pval.pdf',fig_dir),p1,base_height=6)


# Some quality control:
gwas_set=ld_set[snpID==gwas_index,]
rowdata[,snpID:=paste(chr,pos,sep=':')]
rowdata[,snpID:=str_replace(snpID,'chr','')]

# Make heatmap for mvalue at GWAS loci:
gwas_mvalue=mvalue[which(rowdata$snpID%in%gwas_set$snpID),]
pdf(sprintf('%s/mvalue_at_gwas_heatmap.pdf',fig_dir,height=10,width=10))
heatmap.2(as.matrix(gwas_mvalue),Colv=FALSE,Rowv=FALSE,trace='none',dendrogram='none',na.color='black',margin=c(12,12))
dev.off()

# Make boxplot for mvalue at GWAS loci: 
gwas_mvalue_long=melt(gwas_mvalue,variable.name='tissue',value.name='mvalue')
p2=ggplot(gwas_mvalue_long,aes(tissue,mvalue))+geom_boxplot()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
save_plot(sprintf('%s/mvalue_at_gwas_loci.pdf',fig_dir),p2,base_height=8,base_width=8)