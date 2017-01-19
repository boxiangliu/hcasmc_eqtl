# make a list of genes relevant to CAD 
# (and other disease) for running RASQUAL
library(data.table)
library(dplyr)
library(stringr)


# read variants from Nikpay paper: 
nikpay=fread('../processed_data/eCAVIAR/gwas_loci.cad.all.genomewide_fdr_merged.txt')%>%select(chr,pos,markername)


# read variants from latest metabochip paper: 
metabo=fread('../data/gwas/metabochip_novel_lead_variants.txt')%>%select(chr,pos,markername)


# concatenate nikpay and metabochip variants: 
concat=rbind(nikpay,metabo)%>%unique()


# copy position to start and end for foverlap():
concat[,end:=pos]
concat[,start:=pos]
concat$chr=concat$chr%>%str_replace('chr','')%>%paste0('chr',.)


# read genes and their positions: 
gencode=fread('/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf')%>%select(chr=V1,type=V3,start=V4,end=V5,strand=V7,annotation=V9)
gencode=gencode[type=='gene']
gencode[,gene_id:=str_match(annotation,'(?<=gene_id ").+?(?=";)')]
gencode[,gene_name:=str_match(annotation,'(?<=gene_name ").+?(?=";)')]
gencode[,annotation:=NULL]
gencode[,window_start:=ifelse(strand=='+',start-1e6,end-1e6)]
gencode[,window_end:=ifelse(strand=='+',start+1e6,end+1e6)]


# select genes whose TSSes are within 1Mb around any chosen variant:
setkey(concat,chr,start,end)
setkey(gencode,chr,window_start,window_end)
gwas_overlap=foverlaps(concat,gencode)
output=gwas_overlap[,gene_id]%>%unique()


# add some Brian's genes: 
output=c('AHR','CYP1A1','CYP1B1','IL1A','IL1B',output)


# write output:
write.table(output,file="../processed_data/rasqual/prioritized_genes.txt",quote=F,col.names=F,row.names=F)


