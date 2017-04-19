# Library:
library(cowplot)
library(data.table)
library(dplyr)

# Variables: 
out_dir='../processed_data/160715/expression_median_and_rank/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=T)}

# Functions: 
subset_genes=function(x,keep){
	stopifnot(class(keep)=='character')
	y=x[Name%in%keep]
	return(y)
}

median2=function(x,coldata){
	w=data.frame()
	tissue=unique(coldata$tissue)
	for (i in tissue){
		soi=coldata[tissue==i,sample]
		y=x[,colnames(x)%in%soi,with=F]
		z=data.frame(apply(y,1,function(k) {median(k,na.rm=T)}))
		colnames(z)=i
		if (nrow(w)==0){
			w=z
		} else {
			w=cbind(w,z)
		}
	}
	return(w)
}


gid2gname=function(x,rowdata){
	x_id=names(x)
	x_name=rowdata$Description[match(x_id,rowdata$Name)]
	names(x)=x_name
	return(x)
}

list2df=function(x){
	y=list()
	for (i in names(x)){
		y[[i]]=data.frame(gene=names(x[[i]]),esi=x[[i]],tissue=i)
	}
	z=Reduce(rbind,y)
	w=dcast(z,gene~tissue,value.var='esi')
	return(w)
}


# Main: 
# Read RPKM: 
x=fread('../processed_data/160715/combined.gtex.hcasmc.rpkm',header=T)
setnames(x,'9070202_Nextseq','9070202')


# Subset to protein coding gene and lncRNA:
keep=fread('../data/gtex/gencode.v19.genes.v6p.patched_contigs_genetypes.bed',header=F)%>%filter(V1%in%as.character(1:22))%>%filter(V6%in%c('protein_coding','lincRNA'))
x=subset_genes(x,keep$V5)


# Split x into rowdata, coldata and rpkm: 
rowdata=x[,.(Name,Description)]
rpkm=x[,c('Name','Description'):=NULL]
coldata=data.table(sample=colnames(rpkm))


# Read sample name to tissue table: 
sample_list=fread('../data/gtex/SAMPID_SMTSD.sorted.filtered.with_hcasmc.txt',header=F)
setnames(sample_list,c('sample','tissue'))
coldata=merge(coldata,sample_list)
coldata=coldata[match(colnames(rpkm),coldata$sample)]


# Calculate median of each tissue and each gene:
tissue_median=median2(rpkm,coldata)
rownames(tissue_median)=rowdata$Name
tissue_rank=apply(tissue_median,2,rank)
tissue_rank=data.table(tissue_rank)


# Output: 
fwrite(tissue_median,sprintf('%s/median.tsv',out_dir),sep='\t')
fwrite(tissue_rank,sprintf('%s/rank.tsv',out_dir),sep='\t')
