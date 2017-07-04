library(data.table)

bed_fn='../processed_data/atacseq/count/intersect.bed'
out_dir='../processed_data/atacseq/count/'

col=c(1,2,3,4,7,10)
bed=fread(bed_fn,select=c(4,11,col+11),col.names=c('name.a','sample','chr','start','end','name.b','signal','peak'))
bed[,c('left','right'):=list(min(start),max(end)),by='name.a']
peak=unique(bed[,list(chr,left,right)])
fwrite(peak,sprintf('%s/merged_peak.bed',out_dir),col.names=FALSE,sep='\t')
