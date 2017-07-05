library(data.table)

in_fn='../processed_data/atacseq/count/merged_peak.bed'
out_dir='../processed_data/atacseq/rasqual/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

bed=fread(in_fn,col.names=c('chr','start','end'))
bed[,c('source','feature','score','strand','frame','id'):=list('.','gene','.','+','.',paste(chr,start,end,sep='_'))]
bed[,attribute:=sprintf('gene_id "%s"; gene_name "%s";',id,id)]
bed[,id:=NULL]

temp=copy(bed)
temp$feature='exon'

bed=rbind(bed,temp)
setorder(bed,chr,start,end)
setcolorder(bed,c('chr','source','feature','start','end','score','strand','frame','attribute'))

fwrite(bed,sprintf('%s/merged_peak.gtf',out_dir),col.names=FALSE,sep='\t',quote=F)