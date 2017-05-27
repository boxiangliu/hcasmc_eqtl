annotation=fread('/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf')
annotation[,gene_id:=str_extract(V9,'(?<=gene_id ")(.+?)(?=";)')]
annotation[,gene_name:=str_extract(V9,'(?<=gene_name ")(.+?)(?=";)')]
annotation[,type:=str_extract(V9,'(?<=transcript_type ")(.+?)(?=";)')]
annotation[,V4:=V4-1]
annotation=annotation[V3=='gene',]
annotation[,paste0('V',c(2,3,6,8,9)):=NULL]
fwrite(annotation,'/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.bed',sep='\t',col.names=F)