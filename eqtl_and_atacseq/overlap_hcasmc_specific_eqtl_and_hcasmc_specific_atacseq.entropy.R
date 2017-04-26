library(data.table)
library(dplyr)
library(dtplyr)
library(stringr)
library(Hmisc)
library(cowplot)
source('gwas_atacseq_overlap/utils.R')

# Variables 
fig_dir='../figures/eqtl_and_atacseq/overlap_hcasmc_specific_eqtl_and_hcasmc_specific_atacseq/'
if (!dir.exists(fig_dir)) dir.create(fig_dir)


# Functions
parse=function(id){
	parsed=str_split_fixed(id,'_',6)[,1:5]%>%as.data.table()
	setnames(parsed,c('fid','chr','pos','ref','alt'))
	if (!str_detect(parsed$chr[1],'chr')){
		parsed=parsed%>%mutate(chr=paste0('chr',chr))
	}
	return(parsed)
}


get_gene_name_and_id=function(gencode_file){
	gencode=fread(gencode_file)
	gencode=gencode%>%filter(V3=="gene")
	gene_id=str_extract(gencode$V9,'(?<=gene_id ")(ENSG.+?)(?=";)')
	gene_name=str_extract(gencode$V9,'(?<=gene_name ")(.+?)(?=";)')
	stopifnot(length(gene_id)==length(gene_name))
	x=data.table(gene_id=gene_id,gene_name=gene_name)
	return(x)
}


id2name=function(id,gene_name_and_id){
	gene_name=gene_name_and_id$gene_name[match(id,gene_name_and_id$gene_id)]
	return(gene_name)
}


subset2bestQTL=function(x,by,rank){
	setnames(x,c(by,rank),c('fid','logpval'))
	x=x%>%group_by(fid)%>%mutate(is_best=(logpval==max(logpval)))
	x=x%>%filter(is_best==TRUE)
	setnames(x,c('fid','logpval'),c(by,rank))
	x[,is_best:=NULL]
	return(as.data.table(x))
}


append_column=function(x,y,id=c('fid','chr','pos'),col='svalue'){
	setnames(x,id,c('fid','chr','pos'))
	setnames(y,id,c('fid','chr','pos'))
	setnames(y,col,'s')
	if (!is.character(x$pos)) {x[,pos:=as.character(pos)]}
	if (!is.character(y$pos)) {y[,pos:=as.character(pos)]}
	x[,tmp_id:=paste(fid,chr,pos,sep='_')]
	y[,tmp_id:=paste(fid,chr,pos,sep='_')]
	x$new=y[match(x$tmp_id,y$tmp_id),s]
	setnames(x,'new',col)
	x[,tmp_id:=NULL]
	y[,tmp_id:=NULL]
	x[,pos:=as.integer(pos)]
	y[,pos:=as.integer(pos)]
	setnames(x,c('fid','chr','pos'),id)
	setnames(y,c('fid','chr','pos'),id)
	setnames(y,'s',col)
	return(x)
}




## variables
chunk_size=1e6
in_file='../processed_data/eqtl_and_atacseq/specificity.mean.txt'
eqtl_file='../data/eQTL/rasqual/expressedGenes.padj.txt'
fig_dir='../figures/eqtl_and_atacseq/'
gencode_file='/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf'
atacseq_file='../data/atacseq/fbs/2305/out/peak/idr/optimal_set/2305_ppr.IDR0.1.filt.narrowPeak'
roadmap_dir='/mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/'


## main
# Read eQTL:
eqtl=fread(paste0("cat ", eqtl_file, ' | cut -f1,2,3,4,5,6,7,12,25,26,29'))
eqtl=eqtl[r2_rsnp>0.8,]
eqtl[,logpval:=-log10(pval)]


# Subset to bestQTL:
bestQTL=subset2bestQTL(x=eqtl,by='fid',rank='logpval')
bestQTL=bestQTL%>%mutate(start=pos,end=pos)
bestQTL[,id:=paste(fid,chr,pos,sep="_")]


# Get gene id and gene name: 
gene_name_and_id=get_gene_name_and_id(gencode_file)


# Read HCASMC eQTL specificity: 
i=0
container=list()
eof=FALSE
while(!eof){
	# Read current chunk:
	q=fread(in_file,nrows=chunk_size,skip=chunk_size*i)
	setnames(q,c('id','q'))
	q=cbind(parse(q$id),q)
	q$fid=id2name(q$fid,gene_name_and_id)
	q[,id:=paste(fid,chr,pos,sep='_')]
	container[[length(container)+1]]=q[which(id%in%bestQTL$id),]
	i=i+1
	message('chunk ',i)
	eof=ifelse(nrow(q)==chunk_size,FALSE,TRUE)
}


# Concatenate all chunks:
q=Reduce(rbind,container)


# Add HCASMC eQTL specificity (Q) to bestQTL: 
bestQTL=append_column(bestQTL,q,col='q')
bestQTL=bestQTL[!is.na(q)]


# Normalize Q: 
bestQTL[,qnorm:=1-q/max(q,na.rm=T)]


# Bin Q value:
bestQTL[,specific:=cut2(qnorm,g=2)]


# Read peak specificity index:
peak_sp=fread('../processed_data/hcasmc_specific_open_chromatin/raw_peak_specificity/HCASMC.bed')


# Overlap QTL and peak: 
setkey(peak_sp,chr,start,end)
setkey(bestQTL,chr,start,end)
qtl_ol_peak=foverlaps(bestQTL,peak_sp,nomatch=0)


# Plot QSI and PSI: 
p1=ggplot(qtl_ol_peak,aes(psi,qnorm))+geom_point()+stat_smooth(method='lm')+xlab("Peak specificity index")+ylab('QTL specificity index')
save_plot(sprintf("%s/psi_vs_qsi.pdf",fig_dir),p1)


# Fit linear regression on QTL specificity index and peak specificity index:
fit=lm(qnorm~psi,data=qtl_ol_peak)
summary(fit)


# Look at an example (low PSI and low QSI):
qtl_ol_peak[psi<0.5&qnorm<0.25,]
