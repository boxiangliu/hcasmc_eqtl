## library:
library(data.table)
library(dplyr)
library(dtplyr)
library(stringr)
library(cowplot)
source('gwas_atacseq_overlap/utils.R')


## Function: 
subset2bestQTL=function(x,by,rank){
	setnames(x,c(by,rank),c('fid','logpval'))
	x=x%>%group_by(fid)%>%mutate(is_best=(logpval==max(logpval)))
	x=x%>%filter(is_best==TRUE)
	setnames(x,c('fid','logpval'),c(by,rank))
	x[,is_best:=NULL]
	return(as.data.table(x))
}

get_tss=function(gencode_file){
	gencode=fread(gencode_file)
	gencode=gencode%>%filter(V3=="gene")
	x=gencode%>%mutate(chr=V1,tss=ifelse(V7=='+',V4,V5),gene_id=str_extract(V9,'(?<=gene_id ")(ENSG.+?)(?=";)'),gene_name=str_extract(gencode$V9,'(?<=gene_name ")(.+?)(?=";)'))%>%select(chr,tss,gene_id,gene_name)
	return(x)
}


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


calc_dist2tss=function(seq,tss){
	seq[,center:=round((start+end)/2)]
	dist=outer(seq$center,tss$tss,`-`) # takes a while
	same_chr=outer(seq$chr,tss$chr,`==`)
	dist[!same_chr]=Inf
	stopifnot(nrow(dist)==nrow(seq),ncol(dist)==nrow(tss))
	dist2nearestTSS=apply(dist,1,function(x){x[which.min(abs(x))]})
	return(dist2nearestTSS)
}


shuffle=function(seq,tss){
	seq[,center:=round((start+end)/2)]
	seq[,length:=round(abs(start-end)/2)]
	idx=sample(1:nrow(tss),size=nrow(seq),replace=T)
	shuf=tss[idx,]
	shuf$center=shuf$tss+seq$dist
	shuf$start=shuf$center-seq$length
	shuf$end=shuf$center+seq$length
	shuf=shuf%>%select(chr,start,end)
	return(shuf)
}

standardize=function(x,target){
	x[,length:=round(abs(end-start)/2)]
	target[,length:=round(abs(end-start)/2)]
	total_x=sum(x$length)
	total_target=sum(target$length)
	ratio=total_target/total_x
	x[,center:=round((end+start)/2)]
	x[,target_length:=round(length*ratio)]
	x[,target_start:=center-target_length]
	x[,target_end:=center+target_length]
	x=x%>%select(chr,start=target_start,end=target_end)
	return(x)
}


define_peak_region=function(x,size=100){
	x[,peak_pos:=start+peak]
	x[,start:=peak_pos-size]
	x[,end:=peak_pos+size]
	x[,peak_pos:=NULL]
	return(x)
}

run_overlap=function(x,tissue,qtl=bestQTL,base=atacseq,breaks=c(seq(1,6),Inf)){
	y=copy(x)
	tryCatch(setnames(y,c('chrom','chromStart','chromEnd'),c('chr','start','end')),error=function(e){message("Can't rename. Already renamed?")})
	y=define_peak_region(y)

	# Overlap significant eQTL with ATACseq peaks: 
	setkey(qtl,chr,start,end)
	setkey(y,chr,start,end)
	overlapped=foverlaps(qtl,y)

	pct=plot_overlap(overlapped,ratio=nrow(y)/nrow(base),breaks=breaks,plot.it=F)
	pct$tissue=tissue
	return(pct)
}


## variables: 
eqtl_file='../data/eQTL/rasqual/expressedGenes.padj.txt'
atacseq_file='../data/atacseq/fbs/2305/out/peak/idr/optimal_set/2305_ppr.IDR0.1.filt.narrowPeak'
gencode_file='/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf'
roadmap_dir='/mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/'
tmp_dir='/srv/scratch/bliu2/HCASMC_eQTL/eqtl_and_atacseq/'
fig_dir='../figures/eqtl_and_atacseq/'
if (!dir.exists(tmp_dir)){dir.create(tmp_dir,recursive=T)}
if (!dir.exists(fig_dir)){dir.create(fig_dir,recursive=T)}

## main: 
# Read eQTL and ATACseq data: 
eqtl=fread(paste0("cat ", eqtl_file, ' | cut -f1,2,3,4,5,6,7,12,25,26,29'))
atacseq=fread(atacseq_file)


# Keep eQTLs with rsnp r2 > 0.8, subset to best eQTL per gene, and select significant eQTL (FDR < 0.05): 
eqtl=eqtl[r2_rsnp>0.8,]
eqtl=eqtl%>%mutate(start=pos,end=pos)
eqtl[,logpval:=-log10(pval)]
bestQTL=subset2bestQTL(x=eqtl,by='fid',rank='logpval')
sig=eqtl[qval<0.05,]


# Select only unique ATACseq peaks: 
setnames(atacseq,c('chr','start','end','name','score','strand','signal','negLogPval','negLogQval','peak'))
atacseq=atacseq%>%select(chr,start,end)%>%unique()


# Overlap eQTL and ATACseq peaks: 
setkey(bestQTL,chr,start,end)
setkey(atacseq,chr,start,end)
overlapped=foverlaps(bestQTL,atacseq)


# Plot the percentage of eQTL overlapping with ATACseq:
pdf(paste0(fig_dir,'/eqtl_overlap_atacseq.pdf'))
plot_overlap(overlapped)
dev.off()


# Generate random genomic sequences (matched by size and distance to TSS with the ATACseq peaks): 
tss=get_tss(gencode_file)
tss=tss[gene_name%in%unique(eqtl$fid),]
atacseq$dist2tss=calc_dist2tss(atacseq,tss)
rndseq=shuffle(atacseq,tss)
pdf('6.pdf')
hist(atacseq$dist2tss,xlim=c(-1e6,1e6),breaks=10000)
dev.off()

# Overlap eQTL and random genomic sequences: 
setkey(rndseq,chr,start,end)
overlapped_rnd=foverlaps(bestQTL,rndseq)


# Plot the percentage of eQTL overlapping with random genomic sequence:
pdf(paste0(fig_dir,'/eqtl_overlap_random_sequences.pdf'))
plot_overlap(overlapped_rnd)
dev.off()


# Read in Roadmap DNase-seq data: 
atac_gastric=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E094-DNase.macs2.narrowPeak.gz')
# atac_ovary=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E097-DNase.macs2.narrowPeak.gz')
# atac_pancreas=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E098-DNase.macs2.narrowPeak.gz')
# atac_psoasMuscle=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E100-DNase.macs2.narrowPeak.gz')
# atac_smallIntestine=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E109-DNase.macs2.narrowPeak.gz')


# Standardize DNase-seq data to have the same total length as the ATACseq data: 
atac_gastric=atac_gastric%>%select(chr=chrom,start=chromStart,end=chromEnd)%>%unique()
atac_gastric=standardize(atac_gastric,atacseq)
atacseq[,length:=abs(end-start)]
atac_gastric[,length:=abs(end-start)]


# Overlap gastric DNase-seq with bestQTL: 
setkey(atac_gastric,chr,start,end)
overlapped_gastric=foverlaps(bestQTL,atac_gastric)


# Plot the percentage of eQTL overlapping with gastric DNase-seq:
pdf(paste0(fig_dir,'/eqtl_overlap_gastric.pdf'))
plot_overlap(overlapped_gastric)
dev.off()


# Overlap significant eQTL (FDR<0.05) and ATACseq peaks: 
setkey(sig,chr,start,end)
setkey(atacseq,chr,start,end)
overlapped=foverlaps(sig,atacseq)


# Plot the percentage of eQTL overlapping with ATACseq:
pdf(paste0(fig_dir,'/sig_eqtl_overlap_atacseq.pdf'))
plot_overlap(overlapped)
dev.off()


# Overlap gastric DNase-seq with bestQTL: 
setkey(atac_gastric,chr,start,end)
overlapped_gastric=foverlaps(sig,atac_gastric)


# Plot the percentage of eQTL overlapping with gastric DNase-seq:
pdf(paste0(fig_dir,'/sig_eqtl_overlap_gastric.pdf'))
plot_overlap(overlapped_gastric)
dev.off()


# Redefine ATACseq regions: 
atac_hcasmc=fread(atacseq_file)
setnames(atac_hcasmc,c('chr','start','end','name','score','strand','signal','negLogPval','negLogQval','peak'))
atac_hcasmc=define_peak_region(atac_hcasmc)



# Overlap significant eQTL with ATACseq peaks: 
setkey(bestQTL,chr,start,end)
setkey(atacseq,chr,start,end)
overlapped=foverlaps(bestQTL,atacseq)


pct_hcasmc=plot_overlap(overlapped,plot.it=F)
pct_hcasmc$tissue='HCASMC'


# Read Gastric DHS peaks: 
# atac_gastric=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E094-DNase.macs2.narrowPeak.gz')
# atac_ovary=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E097-DNase.macs2.narrowPeak.gz')
# atac_pancreas=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E098-DNase.macs2.narrowPeak.gz')
# atac_psoasMuscle=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E100-DNase.macs2.narrowPeak.gz')
# atac_smallIntestine=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E109-DNase.macs2.narrowPeak.gz')

DNase_samples=list.files(roadmap_dir,pattern='*DNase.macs2.narrowPeak*',full.names=T)

# Read Roadmap DHS samples:
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


# Plot the overlap between eQTL and open chromatin regions: 
pdf(paste0(fig_dir,'eqtl_overlap_open_chromatin.pdf'))
ggplot(pct[tissue!='HCASMC'],aes(x=as.factor(lb),y=pct*100))+geom_boxplot(outlier.shape=NA)+xlab('-log10(P-value)')+ylab('Percent overlap')+geom_point(data=pct[tissue=='HCASMC'],aes(x=as.factor(lb),y=pct*100),color='red')
dev.off()


# Remove temporary files: 
# unlink(tmp_dir,recursive=T)
