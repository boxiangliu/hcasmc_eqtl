library(data.table)
library(stringr)
args=commandArgs(T)
in_fn=args[1]
out_fn=args[2]
x=fread(in_fn,select=c(1,3,5,6,7),col.names=c('chr','SNP_pos','intron_start','intron_end','intron'))
x[,dist:=pmin(abs(SNP_pos-intron_start),abs(SNP_pos-intron_end))]
x[,dist_nearest_intron:=min(dist),by=c('chr','SNP_pos')]

y=unique(x[dist==dist_nearest_intron,list(chr,SNP_pos,dist_nearest_intron)])
y[,snpID:=paste0(str_replace(chr,'chr',''),":",SNP_pos)]

fwrite(y[,list(snpID,dist_nearest_intron)],out_fn,sep='\t')