# make a list of genes relevant to CAD 
# (and other disease) for running RASQUAL


# library:
library(dplyr)
library(data.table)
library(jsonlite)
library(stringr)


# read GWAS catelog: 
catalog=fread('/srv/persistent/bliu2/HCASMC_eQTL/data/gwas/gwas_catalog_v1.0-associations_e86_r2016-10-23.tsv')


# how many unique SNPs:
length(unique(catalog[,SNPS])) # 24986


# how many unique SNPs with p < 5e-8: 
length(unique(catalog[`P-VALUE`<1e-8,SNPS])) # 6985


# how many traits: 
length(table(catalog[,`DISEASE/TRAIT`])) # 1601


# trait I care about:
# diabetes
# aneurysms
# blood pressure
# cardiac arrest 
# Stroke
# ...


# compare current and past rs id: 
# rsid=catalog[,.(SNPS,SNP_ID_CURRENT)]%>%mutate(SNP_ID_CURRENT=paste0('rs',SNP_ID_CURRENT))
# SNP_ID_CURRENT and SNPS maps to the same coordinate using ENSEMBL. 
# can check with the following command, just change the rsid: 
# response=fromJSON(system("curl 'http://grch37.rest.ensembl.org/variation/human/rs1047891?' -H 'Content-type:application/json'",intern=T))$mappings%>%select(chr=seq_region_name,start,end)

# there are some SNPs that has ';':
# catalog[which(str_detect(rsid$SNPS,';'))[1],]
# need to remove these SNPs. 


# subset the catalog: 
catalog_bak=catalog
catalog=catalog_bak[,.(`DISEASE/TRAIT`,CHR_ID,CHR_POS,`STRONGEST SNP-RISK ALLELE`,SNPS,SNP_ID_CURRENT,`P-VALUE`,PVALUE_MLOG,`OR or BETA`,CNV)]


# there are rows with ambiguous rs IDs: 
# table(nchar(catalog$SNPS))
# catalog[nchar(catalog$SNPS)==0]
# catalog[nchar(catalog$SNPS)==2]
# catalog[nchar(catalog$SNPS)==4]
# catalog[nchar(catalog$SNPS)==6]
# catalog[nchar(catalog$SNPS)==21]


# select SNPs with unambiguous rs ID:
catalog_filt=catalog[which(!is.na(str_match(catalog$SNPS,'^rs[0-9]+$'))),]

# get hg19 coordinates for SNPs: 
collector=data.frame()
for (i in 1:nrow(catalog_filt)){
	message(i)
	rsid=catalog_filt$SNPS[i]
	response=tryCatch({
		tmp=fromJSON(system(sprintf("curl 'http://grch37.rest.ensembl.org/variation/human/%s?' -H 'Content-type:application/json'",rsid),intern=T));
		response=tmp$mappings%>%select(chr=seq_region_name,start,end);
		response$rsid=tmp$name;
		response},
		error=function(e) {
			message(paste(rsid,'cannot be found in ENSEMBL')) 
			return(data.frame(chr=NA,start=NaN,end=NaN))})
	collector=rbind(collector,response)
}


# count how many rows with chromosome names 1 to 22 and X: 
# sum(table(catalog$CHR_ID)[which(names(table(catalog$CHR_ID)) %in% c(as.character(seq(1,22)),'X'))]) # 31997


# select SNPs with chromosome names 1 to 22 and X: 
catalog_filt=catalog[which(catalog$CHR_ID%in%c(as.character(seq(1,22)),'X')),]


# output bed file:
bed=catalog_filt%>%mutate(start=as.integer(CHR_POS)-1,chr=paste0('chr',CHR_ID))%>%select(chr,start,end=CHR_POS)
write.table(bed,file='../processed_data/rasqual/tmp/gwas_catalog.hg38.bed',quote=F,sep="\t",row.names=F,col.names=F)


# liftOver from hg38 to hg19: 
system('/srv/persistent/bliu2/tools/ucsc_tools/liftOver ../processed_data/rasqual/tmp/gwas_catalog.hg38.bed /srv/persistent/bliu2/shared/chain_files/hg38ToHg19.over.chain.gz ../processed_data/rasqual/tmp/gwas_catalog.hg19.bed ../processed_data/rasqual/tmp/gwas_catalog.unmapped.bed')


# read liftOver result:
catalog_hg19=fread('../processed_data/rasqual/tmp/gwas_catalog.hg19.bed')%>%select(chr=V1,pos=V3)
catalog_hg19$study='nhgri'

# read Nikpay paper: 
nikpay=fread('../processed_data/eCAVIAR/gwas_loci.cad.all.genomewide_fdr_merged.txt')%>%select(chr,pos)
nikpay$study='nikpay'


# merge GWAS catalog and Nikpay: 
merged=unique(rbind(nikpay,catalog_hg19))
merged[,end:=pos]
merged[,start:=pos]
merged$chr=merged$chr%>%str_replace('chr','')%>%paste0('chr',.)


# read genes and their positions: 
gencode=fread('/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf')%>%select(chr=V1,type=V3,start=V4,end=V5,strand=V7,annotation=V9)
gencode=gencode[type=='gene']
gencode[,gene_id:=str_match(annotation,'(?<=gene_id ").+?(?=";)')]
gencode[,gene_name:=str_match(annotation,'(?<=gene_name ").+?(?=";)')]
gencode[,annotation:=NULL]


# select genes whose TSSes are within 1Mb around any chosen variant:
gencode[,window_start:=ifelse(strand=='+',start-1e6,end-1e6)]
gencode[,window_end:=ifelse(strand=='+',start+1e6,end+1e6)]
setkey(merged,chr,start,end)
setkey(gencode,chr,window_start,window_end)
gwas_overlap=foverlaps(merged,gencode)
output=gwas_overlap[,.(gene_name,study)]%>%unique()%>%arrange(study)
table(output$study)



nikpay[,end:=pos]
nikpay[,start:=pos]
nikpay$chr=nikpay$chr%>%str_replace('chr','')%>%paste0('chr',.)

setkey(nikpay,chr,start,end)
setkey(gencode,chr,window_start,window_end)
gwas_overlap=foverlaps(nikpay,gencode)
output=gwas_overlap[,gene_id]%>%unique()


# add some Brian's genes: 
output=c('AHR','CYP1A1','CYP1B1','IL1A','IL1B',output)


# write output: 
write.table(output,file="../processed_data/rasqual/prioritized_genes.txt",quote=F,col.names=F,row.names=F)


