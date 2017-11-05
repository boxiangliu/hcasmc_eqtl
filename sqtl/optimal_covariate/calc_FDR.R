library(data.table)
library(foreach)
library(dplyr)
library(dtplyr)

args=commandArgs(T)
in_fn=args[1]
out_fn=args[2]
geno_pc=as.integer(args[3])
splice_pc=as.integer(args[4])

message('INFO - reading p-values')
fastqtl=fread(sprintf('zcat %s',in_fn),select=c(1,2,4),col.names=c('cluster','snp','p'))

fastqtl=fastqtl%>%
group_by(cluster)%>%
mutate(pmin=min(p))%>%
filter(p==pmin)

stopifnot(all(fastqtl$p==fastqtl$pmin))


message('INFO - calculating FDR')
padj=p.adjust(fastqtl$p,method='fdr')

message('INFO - counting significant sQTL')
sig=foreach(FDR=c(1e-8,1e-6,1e-4,0.01,0.05),.combine='rbind')%do%{
	data.table(FDR=FDR,sig=sum(padj<FDR,na.rm=TRUE))
}

sig$geno_pc=geno_pc
sig$splice_pc=splice_pc


message('INFO - writing output')
fwrite(sig,out_fn,sep='\t',col.names=FALSE)


message('INFO - done')