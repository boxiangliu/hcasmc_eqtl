## library:
library(data.table)
library(dplyr)
library(dtplyr)
library(stringr)
library(cowplot)

## Function: 
subset2bestQTL=function(x,by,rank){
	setnames(x,c(by,rank),c('fid','logpval'))
	x=x%>%group_by(fid)%>%mutate(is_best=(logpval==max(logpval)))
	x=x%>%filter(is_best==TRUE)
	setnames(x,c('fid','logpval'),c(by,rank))
	x[,is_best:=NULL]
	return(as.data.table(x))
}


## variables: 
eqtl_file='../data/eQTL/rasqual/expressedGenes.padj.txt'
atacseq_file='../data/atacseq/fbs/2305/out/peak/idr/optimal_set/2305_ppr.IDR0.1.filt.narrowPeak'

## main: 
# Read eQTL and ATACseq data: 
eqtl=fread(paste0("cat ", eqtl_file, ' | cut -f1,2,3,4,5,6,7,12,25,26,29'))
atacseq=fread(atacseq_file)

# Keep eQTLs with rsnp r2 > 0.8: 
eqtl=eqtl[r2_rsnp>0.8,]

# Select only unique ATACseq peaks: 
setnames(atacseq,c('chr','start','end','name','score','strand','signal','negLogPval','negLogQval','peak'))
atacseq=atacseq%>%select(chr,start,end)%>%unique()


# Overlap eQTL and ATACseq peaks: 
eqtl=eqtl%>%mutate(start=pos,end=pos)
setkey(eqtl,chr,start,end)
setkey(atacseq,chr,start,end)
overlapped=foverlaps(eqtl,atacseq)
overlapped[,logpval:=-log10(pval)]



# Subset to significant eQTLs (p<0.01):
# sig=overlapped%>%filter(logpval>2)


# If an eQTL regulates multiple genes, count it only once:
# sig_once=sig%>%group_by(chr,start,end,rsid,i.start,i.end)%>%dplyr::summarise(logpval=max(logpval))
# sig_once=as.data.table(sig_once)


# Take a single eQTL per peak to account for LD: 
# na=sig_once%>%filter(is.na(start))
# notNa=sig_once%>%filter(!is.na(start))
# notNa=notNa%>%group_by(chr,start,end)%>%dplyr::summarise(logpval=max(logpval))
# na=na%>%select(chr,start,end,logpval)
# sig_LD=rbind(na,notNa)
# sig_LD=as.data.table(sig_LD)


# Plot the percentage of eQTL overlapping with ATACseq: 
# sig_LD[,bin:=cut2(logpval,g=20)]
# sig_LD[,pct:=mean(!is.na(start)),by='bin']
# pct=sig_LD%>%select(bin,pct)%>%unique()
# pct$lb=as.numeric(str_extract(pct$bin,'(?<=[\\[])(.+?)(?=,)'))
# pct=pct%>%arrange(lb)
# p2=ggplot(pct,aes(x=lb,y=100*pct))+geom_point(size=5,alpha=0.8)+stat_smooth()+xlab('-log10(P-value)')+ylab('Percentage overlap')
# pdf()
# p2
# dev.off()


# Subset to one significant eQTL per gene: 
bestQTL=subset2bestQTL(x=overlapped,by='fid',rank='logpval')
bestQTLbak=bestQTL


# Plot the percentage of eQTL overlapping with ATACseq:
bestQTL[,bin:=cut2(logpval,g=20)]
bestQTL[,pct:=mean(!is.na(start)),by='bin']
pct=bestQTL%>%select(bin,pct)%>%unique()
pct$lb=as.numeric(str_extract(pct$bin,'(?<=[\\[])(.+?)(?=,)'))
pct=pct%>%arrange(lb)
p3=ggplot(pct,aes(x=lb,y=100*pct))+geom_point(size=5,alpha=0.8)+stat_smooth()+xlab('-log10(P-value)')+ylab('Percentage overlap')
pdf()
p3
dev.off()