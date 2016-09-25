# add rs id from dbsnp to eqtl file: 

library(R.utils)
library(data.table)
library(dplyr)
library(stringr)

# command args:
args=commandArgs(T,T)
dbsnp_file=args$dbsnp
eqtl_file=args$eqtl
out_file=args$out
# dbsnp_file='../processed_data/160811/tested_snps_at_egenes.dbsnp146.txt'
# eqtl_file='../processed_data/160811/tested_snps_at_egenes.txt'
# out_file='../processed_data/160811/tested_snps_at_egenes.rsid.txt'

# read dbsnp:
dbsnp=fread(dbsnp_file)
setnames(dbsnp,c('chr','pos','rsid','strand','geno'))


# read eqtl: 
eqtl=fread(eqtl_file)%>%dplyr::select(1:6)
setnames(eqtl,c('pheno','geno','dist','pval','beta','se'))


# add ID column to dbsnp:
dbsnp[,ID:=paste(chr,pos,sep="_")]


# add ID column to eqtl: 
temp=eqtl[,str_split_fixed(geno,pattern='_',n=4)[,c(1,2)]]
eqtl$ID=apply(temp,1,paste,collapse="_")


# add rsid to eqtl: 
idx=match(eqtl$ID,dbsnp$ID)
eqtl$rsid=dbsnp$rsid[idx]
eqtl[,rsid:=ifelse(is.na(rsid),geno,rsid)]


# write output: 
output=eqtl%>%select(-ID)
write.table(output,out_file,quote=F,sep='\t',row.names=F,col.names=F)