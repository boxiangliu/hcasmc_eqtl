## library:
library(data.table)
library(dplyr)
library(dtplyr)
library(stringr)
library(cowplot)
# source('gwas_atacseq_overlap/utils.R')


## Function:
plot_overlap=function(overlapped,ratio=1,g=20,plot.it=T,breaks=c(1,2,3,4,5,6,Inf)){
	# overlapped[,bin:=cut2(logpval,g=g)]
	overlapped[,bin:=cut(logpval,breaks=breaks)]
	overlapped[,pct:=mean(!is.na(start))/ratio,by='bin']
	pct=overlapped%>%select(bin,pct)%>%unique()
	pct$lb=as.numeric(str_extract(pct$bin,'(?<=[\\(\\[])(.+?)(?=,)'))
	pct=pct%>%arrange(lb)
	if (plot.it){
		p=ggplot(pct,aes(x=lb,y=100*pct))+geom_point(size=5,alpha=0.8)+stat_smooth()+xlab('-log10(P-value)')+ylab('Percentage overlap')
		return(p)
	} else {
		return(pct)
	}
}

run_overlap=function(x,tissue,qtl=bestQTL,base=atacseq,breaks=c(seq(1,8),Inf)){
	# Overlap significant eQTL with ATACseq peaks: 
	setkey(qtl,chr,start,end)
	setkey(x,chr,start,end)


	pct=data.frame()
	for (b in breaks){
		y=qtl[logpval>b,]
		overlapped=foverlaps(y,x,nomatch=0)
		n_overlap=nrow(unique(overlapped[,list(fid,rsid,chr,pos,ref,alt)]))
		pct=rbind(pct,data.frame(lb=b,pct=n_overlap/nrow(y)))
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

## variables: 
eqtl_file='../processed_data/compare_rasqual_and_fastqtl/top_snp_per_gene/rasqual.txt'
dhs_dir='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_adult_filt/'
tmp_dir='/srv/scratch/bliu2/HCASMC_eQTL/eqtl_and_atacseq/'
fig_dir='../figures/eqtl_and_atacseq/overlap_eqtl_and_atacseq/'
if (!dir.exists(tmp_dir)){dir.create(tmp_dir,recursive=T)}
if (!dir.exists(fig_dir)){dir.create(fig_dir,recursive=T)}


## main: 
# Read eQTL and ATACseq data:
ras=fread(eqtl_file)
ras[,logpval:=-log10(pval)]
ras[,c('start','end'):=pos]


# Read Roadmap DHS samples:
DNase_samples=list.files(dhs_dir,pattern='bed')
atacseq=list()
for (i in 1:length(DNase_samples)){
	in_fn=DNase_samples[i]
	sample=str_replace(in_fn,'.merged.bed','')
	print(sprintf('INFO - %s',sample))
	tmp=fread(sprintf('%s/%s',dhs_dir,in_fn))
	tmp=tmp[!chr%in%paste0('chr',c('X','Y','M')),]
	tmp[,c('start','end'):=list(start-500,end+500)]
	atacseq[[sample]]=tmp
}



# Calculate percentage overlap between eQTL and open chromatin (without correcting for number of peaks: 
pct=data.frame()
for (i in names(atacseq)){
	pct=rbind(pct,run_overlap(atacseq[[i]],i,qtl=ras,breaks=c(seq(1,8),Inf),base='none'))
}
setDT(pct)
pct=pct[!is.na(pct),]

# Plot the overlap between eQTL and open chromatin regions: 
pdf(sprintf('%s/overlap_eqtl_and_atacseq.pdf',fig_dir))
ggplot(pct[tissue!='HCASMC'],aes(x=as.factor(lb),y=pct*100))+geom_boxplot(outlier.shape=NA)+xlab('-log10(P-value)')+ylab('Percent overlap')+geom_line(data=pct[tissue=='HCASMC'],aes(x=lb,y=pct*100),color='red')+geom_point(data=pct[tissue=='HCASMC'],aes(x=lb,y=pct*100),size=5,color='red')
dev.off()

# Downsample peaks:
set.seed(42)
atacseq_ds=downsample(atacseq)

# Calculate percentage overlap between eQTL and open chromatin (without correcting for number of peaks: 
pct=data.frame()
for (i in names(atacseq_ds)){
	pct=rbind(pct,run_overlap(atacseq_ds[[i]],i,qtl=ras,breaks=c(seq(1,8),Inf),base='none'))
}
setDT(pct)
pct=pct[!is.na(pct),]

# Plot the overlap between eQTL and open chromatin regions: 
pdf(sprintf('%s/overlap_eqtl_and_atacseq.downsample.pdf',fig_dir))
ggplot(pct[tissue!='HCASMC'],aes(x=as.factor(lb),y=pct*100))+geom_boxplot(outlier.shape=NA)+xlab('-log10(P-value)')+ylab('Percent overlap')+geom_line(data=pct[tissue=='HCASMC'],aes(x=lb,y=pct*100),color='red')+geom_point(data=pct[tissue=='HCASMC'],aes(x=lb,y=pct*100),size=5,color='red')
dev.off()


# Remove temporary files: 
unlink(tmp_dir,recursive=T)
