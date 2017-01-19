library(data.table)
library(stringr)
library(dplyr)
library(dtplyr)

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

check_variant_alleles=function(chr,pos,ref,alt,vcf_file){
	correct_alleles=data.frame()
	for (i in 1:length(ref)){
		record=system(sprintf('tabix %s %s:%s-%s',vcf_file,chr[i],pos[i],pos[i]),intern=T)
		record=str_split(record,'\t')[[1]][1:5]
		if (record[4]!=ref[i] || record[5]!=alt[i]){message(sprintf("[%s:%s;%s] (%s,%s) -> (%s,%s)",chr[i],pos[i],record[3],ref[i],alt[i],record[4],record[5]))}
		correct_alleles=rbind(correct_alleles,data.frame(ref=record[4],alt=record[5]))
	}
	return(correct_alleles)
}


# variables:
in_file='../processed_data/mpra/naive_overlap_eQTL/gwas_eQTL_naive_overlap.txt'
out_file='../processed_data/mpra/forNathan/naive_overlap_variants.bed'
gencode_file='/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf'
vcf_file='../processed_data/160604_phasing/phased_and_imputed_gprobs/phased_and_imputed.all.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz'


# main:
variants=fread(in_file)
gene_name_and_id=get_gene_name_and_id(gencode_file)
variants$gene_id=name2id(variants$fid,gene_name_and_id)
correct_alleles=check_variant_alleles(variants$chr,variants$pos,variants$ref,variants$alt,vcf_file)
variants$ref=correct_alleles$ref
variants$alt=correct_alleles$alt
output=variants%>%filter(rsq_rsnp>0.8)%>%mutate(start=pos-1)%>%select(chr,start,end=pos,ref,alt,gene_id,rsid,pval,padj,rank)
if(!dir.exists(dirname(out_file))){dir.create(dirname(out_file))}
write.table(output,out_file,sep='\t',col.names=T,row.names=F,quote=F)



