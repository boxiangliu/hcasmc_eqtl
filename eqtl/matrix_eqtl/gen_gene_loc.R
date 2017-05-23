#!/usr/bin/env Rscript
# bosh liu
# 2016/05/08
# durga
# generate gene location file

library(data.table)
library(stringr)

# command args: 
args=commandArgs(T)
gencode_file=args[1]
output_file=args[2]
if (!dir.exists(dirname(output_file))) {dir.create(dirname(output_file),recursive=TRUE)}

# read gene code annotation:
gencode=fread(gencode_file,header=F)


# setnames:
setnames(gencode,c('V1','V2','V3','V4','V5','V9'),c('chr','source','type','left','right','annotation'))



# subset to only genes:
gencode=gencode[type=='gene',]


# remove "chr" from chromosome names: 
# gencode[,chr:=str_replace(chr,"chr","")]


# extract gene_id:
gene_id=str_match(gencode[,annotation],'gene_id \\"(ENSG.+?)\\";')[,2]


# construct gene location data.frame:
gene_loc=data.frame(id=gene_id, chr=gencode[,chr],left=gencode[,left],right=gencode[,right])



# write output:
write.table(gene_loc,file=output_file,row.names=F,quote=F,sep='\t')
