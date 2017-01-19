# library: 
library(data.table)
library(dplyr)
library(stringr)


# functions: 
get_gene_name_and_id=function(gencode_file){
	gencode=fread(gencode_file)
	gencode=gencode%>%filter(V3=="gene")
	gene_id=str_extract(gencode$V9,'(?<=gene_id ")(ENSG.+?)(?=";)')
	gene_name=str_extract(gencode$V9,'(?<=gene_name ")(.+?)(?=";)')
	stopifnot(length(gene_id)==length(gene_name))
	x=data.table(gene_id=gene_id,gene_name=gene_name)
	return(x)
}

name2id=function(name,gene_name_and_id){
	gene_id=c()
	for (n in name){
		x=gene_name_and_id%>%filter(gene_name==n)%>%select(gene_id)%>%unlist()
		if (length(x)>1) {stop("Name is not unique!")}
		gene_id=c(gene_id,x)
	}
	names(gene_id)=NULL
	return(gene_id)
}


# variables:
args=commandArgs(T)
in_file=args[1]
out_file=args[2]
gencode_file='/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf'


# main:
variants=fread(in_file)
gene_name_and_id=get_gene_name_and_id(gencode_file)
output=merge(variants,gene_name_and_id,by='gene_id')
output[,gene_id:=NULL]
write.table(output,out_file,quote=FALSE,row.names=F)