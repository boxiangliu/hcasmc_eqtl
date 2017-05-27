## library
library(data.table)
library(dplyr)
library(dtplyr)
library(stringr)
library(Hmisc)
library(cowplot)
source('gwas_atacseq_overlap/utils.R')

## function
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


define_peak_region=function(x,size=100){
	x[,peak_pos:=start+peak]
	x[,start:=as.integer(peak_pos-size)]
	x[,end:=as.integer(peak_pos+size)]
	x[,peak_pos:=NULL]
	return(x)
}


plot_overlap=function(overlapped,ratio=1,g=20,plot.it=T,breaks=c(0,1,2,3,4,5,6,Inf)){
	# overlapped[,bin:=cut2(logpval,g=g)]
	overlapped[,bin:=cut(logpval,breaks=breaks)]
	overlapped[,pct:=mean(!is.na(start))/ratio,by=c('bin','specific')]
	pct=overlapped%>%select(bin,specific,pct)%>%unique()
	pct$lb=as.numeric(str_extract(pct$bin,'(?<=[\\(\\[])(.+?)(?=,)'))
	pct=pct%>%arrange(lb)
	if (plot.it){
		p=ggplot(pct,aes(x=lb,y=100*pct))+geom_point(size=5,alpha=0.8)+stat_smooth()+xlab('-log10(P-value)')+ylab('Percentage overlap')
		return(p)
	} else {
		return(pct)
	}
}


run_overlap=function(x,tissue,qtl=bestQTL,base=atacseq,breaks=c(seq(1,6),Inf)){
	y=copy(x)
	tryCatch(setnames(y,c('chrom','chromStart','chromEnd'),c('chr','start','end')),error=function(e){message("Can't rename. Already renamed?")})
	y=define_peak_region(y)

	# Overlap significant eQTL with ATACseq peaks: 
	setkey(qtl,chr,start,end)
	setkey(y,chr,start,end)
	overlapped=foverlaps(qtl,y)

	if (base=='none'){
		pct=plot_overlap(overlapped,ratio=1,breaks=breaks,plot.it=F)
	} else {
		pct=plot_overlap(overlapped,ratio=nrow(y)/nrow(base),breaks=breaks,plot.it=F)
	}
	pct$tissue=tissue
	return(pct)
}


downsample=function(x){
	min_num_peaks=min(sapply(x,nrow))
	downsample=list()
	for (i in names(x)){
		downsample[[i]]=x[[i]][sample(1:nrow(x[[i]]),min_num_peaks),]
	}
	return(downsample)
}


count_num_eQTL=function(pct,qtl=bestQTL,breaks=c(0,1,2,3,4,5,6,Inf)){
	qtl1=copy(qtl)
	qtl1[,bin:=cut(logpval,breaks=breaks)]
	label=qtl1%>%group_by(bin,specific)%>%summarise(num_eQTL=n())
	label$tissue='HCASMC'
	pct1=merge(pct,label,by=c('bin','specific','tissue'),all.x=T)
	pct1[,num_eQTL:=ifelse(is.na(num_eQTL),'',num_eQTL)]
	return(pct1)
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
nrow(q)

# Add HCASMC eQTL specificity (Q) to bestQTL: 
bestQTL=append_column(bestQTL,q,col='q')
bestQTL=bestQTL[!is.na(q)]


# Normalize Q: 
bestQTL[,qnorm:=1-q/max(q,na.rm=T)]


# Bin Q value:
bestQTL[,specific:=cut2(qnorm,g=2)]


# Plot HCASMC specific score (Q):
p=ggplot(bestQTL,aes(qnorm))+geom_histogram(binwidth=0.02)+geom_vline(aes(xintercept=median(qnorm)),color='red',linetype=2)+xlab('eQTL specificity')+annotate(geom='text',label=c('HCASMC-specific','Shared'),x=c(0.8,0.1),y=c(3000,3000))
save_plot(paste(fig_dir,'eqtl_specificity_score.pdf'),p)


# Read ATACseq data:
atac_hcasmc=fread(atacseq_file)%>%setnames(.,c('chr','start','end','name','score','strand','signalValue','pValue','qValue','peak'))


# Read Roadmap DHS samples:
DNase_samples=list.files(roadmap_dir,pattern='*DNase.macs2.narrowPeak*',full.names=T)
atacseq=list(atac_hcasmc)
for (i in 1:length(DNase_samples)){
	atacseq[[length(atacseq)+1]]=readAtac(sprintf('zcat %s',DNase_samples[i]))
}


# Assign tissue names:
tissues=c('HCASMC',str_match(DNase_samples,'E[0-9]{3}'))
names(atacseq)=tissues


# Calculate percentage overlap between eQTL and open chromatin: 
pct=data.frame()
for (i in names(atacseq)){
	pct=rbind(pct,run_overlap(atacseq[[i]],i,base=atacseq[['HCASMC']]))
}
pct=pct[!is.na(bin),]
pct=count_num_eQTL(pct)
pct[,specific:=ifelse(as.numeric(specific)==2,'HCASMC-specific','Shared')]


# Plot overlap between eQTL and open chromatin, stratified by HCASMC-specificity: 
p1=ggplot(pct[tissue!='HCASMC'],aes(x=as.factor(lb),y=pct*100))+geom_boxplot(outlier.shape=NA)+xlab('-log10(P-value)')+ylab('Percent overlap')+geom_line(data=pct[tissue=='HCASMC'],aes(x=lb,y=pct*100),color='red')+geom_point(data=pct[tissue=='HCASMC'],aes(x=lb,y=pct*100),size=5,color='red')+geom_text(data=pct[tissue=='HCASMC'],aes(x=lb,y=pct*100,label=num_eQTL))+facet_grid(.~specific)+theme_bw()
save_plot(paste0(fig_dir,'eqtl_overlap_open_chromatin_by_specificity.2quantiles.pdf'),p1)


# Bin eQTL into three groups: 
bestQTL[,specific:=cut2(qnorm,g=3)]

# Calculate percentage overlap between eQTL and open chromatin: 
pct=data.frame()
for (i in names(atacseq)){
	pct=rbind(pct,run_overlap(atacseq[[i]],i,base=atacseq[['HCASMC']]))
}
pct=pct[!is.na(bin),]
pct=count_num_eQTL(pct)
pct=pct[as.numeric(specific)!=2,]
pct[,specific:=ifelse(as.numeric(specific)==3,'HCASMC-specific','Shared')]


# Plot overlap between eQTL and open chromatin, stratified by HCASMC-specificity: 
p2=ggplot(pct[tissue!='HCASMC'],aes(x=as.factor(lb),y=pct*100))+geom_boxplot(outlier.shape=NA)+xlab('-log10(P-value)')+ylab('Percent overlap')+geom_line(data=pct[tissue=='HCASMC'],aes(x=lb,y=pct*100),color='red')+geom_point(data=pct[tissue=='HCASMC'],aes(x=lb,y=pct*100),size=5,color='red')+geom_text(data=pct[tissue=='HCASMC'],aes(x=lb,y=pct*100,label=num_eQTL))+facet_grid(.~specific)+theme_bw()
save_plot(paste0(fig_dir,'eqtl_overlap_open_chromatin_by_specificity.3quantiles.pdf'),p2)



# Downsample:
set.seed(42)
atacseq_ds=downsample(atacseq)


# Bin Q value:
bestQTL[,specific:=cut2(qnorm,g=2)]


# Calculate percentage overlap between eQTL and open chromatin: 
pct=data.frame()
for (i in names(atacseq_ds)){
	pct=rbind(pct,run_overlap(atacseq_ds[[i]],i,base='none'))
}
pct=pct[!is.na(bin),]
pct=count_num_eQTL(pct)
pct[,specific:=ifelse(as.numeric(specific)==2,'HCASMC-specific','Shared')]


# Plot overlap between eQTL and open chromatin, stratified by HCASMC-specificity: 
p3=ggplot(pct[tissue!='HCASMC'],aes(x=as.factor(lb),y=pct*100))+geom_boxplot(outlier.shape=NA)+xlab('-log10(P-value)')+ylab('Percent overlap')+geom_line(data=pct[tissue=='HCASMC'],aes(x=lb,y=pct*100),color='red')+geom_point(data=pct[tissue=='HCASMC'],aes(x=lb,y=pct*100),size=5,color='red')+geom_text(data=pct[tissue=='HCASMC'],aes(x=lb,y=pct*100,label=num_eQTL))+facet_grid(.~specific)+theme_bw()
save_plot(paste0(fig_dir,'eqtl_overlap_open_chromatin_by_specificity.downsample.2quantiles.pdf'),p3)


# Bin Q value:
bestQTL[,specific:=cut2(qnorm,g=3)]


# Calculate percentage overlap between eQTL and open chromatin: 
pct=data.frame()
for (i in names(atacseq_ds)){
	pct=rbind(pct,run_overlap(atacseq_ds[[i]],i,base='none'))
}
pct=pct[!is.na(bin),]
pct=count_num_eQTL(pct)
pct=pct[as.numeric(specific)!=2,]
pct[,specific:=ifelse(as.numeric(specific)==3,'HCASMC-specific','Shared')]


# Plot overlap between eQTL and open chromatin, stratified by HCASMC-specificity: 
p4=ggplot(pct[tissue!='HCASMC'],aes(x=as.factor(lb),y=pct*100))+geom_boxplot(outlier.shape=NA)+xlab('-log10(P-value)')+ylab('Percent overlap')+geom_line(data=pct[tissue=='HCASMC'],aes(x=lb,y=pct*100),color='red')+geom_point(data=pct[tissue=='HCASMC'],aes(x=lb,y=pct*100),size=5,color='red')+geom_text(data=pct[tissue=='HCASMC'],aes(x=lb,y=pct*100,label=num_eQTL))+facet_grid(.~specific)+theme_bw()
save_plot(paste0(fig_dir,'eqtl_overlap_open_chromatin_by_specificity.downsample.3quantiles.pdf'),p4)