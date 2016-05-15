#!/usr/bin/env Rscript
# bosh liu
# 2016/05/08
# durga
# generate gene location file


# paths: 
gencode_file='/srv/persistent/bliu2/shared/annotation/gencode/v14/gencode.v14.annotation.gtf'
output_dir='../processed_data/031_gen_gene_loc'


# read gene code annotation:
gencode=fread(gencode_file,header=F)


# setnames:
setnames(gencode,c('V1','V2','V3','V4','V5','V9'),c('chr','source','type','left','right','annotation'))



# subset to only genes:
gencode=gencode[type=='gene',]


# remove "chr" from chromosome names: 
gencode[,chr:=str_replace(chr,"chr","")]


# extract gene_id:
gene_id=str_match(gencode[,annotation],'gene_id \\"(ENSG.+?)\\";')[,2]


# construct gene location data.frame:
gene_loc=data.frame(id=gene_id, chr=gencode[,chr],left=gencode[,left],right=gencode[,right])



# write output:
write.table(gene_loc,file=sprintf('%s/gene_loc.txt',output_dir),row.names=F,quote=F,sep='\t')
