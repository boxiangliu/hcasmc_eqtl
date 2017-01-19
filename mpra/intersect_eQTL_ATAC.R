library(data.table)
library(dplyr)
library(dtplyr)
library(cowplot)

eqtl_file='../processed_data/mpra/naive_overlap_eQTL/gwas_eQTL_naive_overlap.txt'
atacseq_file='../data/atacseq/fbs/2305/out/peak/idr/optimal_set/2305_ppr.IDR0.1.filt.narrowPeak'

eqtl=fread(eqtl_file)
unique(eqtl[,.(rsid,chr,pos)])
eqtl
atacseq=fread(atacseq_file)
setnames(atacseq,c('chr','start','end','name','score','strand','signalValue','pValue','qValue','peak'))

eqtl=eqtl%>%mutate(start=pos)
eqtl=eqtl%>%mutate(end=pos)
eqtl=eqtl[rsq_rsnp>=0.8,]

setkey(eqtl,chr,start,end)
setkey(atacseq,chr,start,end)
overlap=foverlaps(eqtl,atacseq)
overlap[!is.na(peak),.(fid,rsid,pval,padj,rank)]%>%as.data.frame()%>%arrange(rank)%>%unique()
overlap[!is.na(peak),.(fid,rsid,pval,padj,rank)]%>%as.data.frame()%>%arrange(rank)%>%unique()%>%select(padj)%>%unlist()%>%median()
overlap[is.na(peak),.(fid,rsid,pval,padj,rank)]%>%as.data.frame()%>%arrange(rank)%>%unique()%>%filter(padj<0.0034931)%>%nrow()
eqtl%>%arrange(rank)%>%as.data.frame()%>%head(10)
ggplot(overlap,aes(y=-log10(pval),x=is.na(peak)))+geom_violin()+geom_boxplot(width=0.1)

