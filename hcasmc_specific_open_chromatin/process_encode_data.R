# Select ENCODE cell lines with no "AUDIT NOT COMPLIANT" or "AUDIT ERROR" flag.
# Subsample each selected sample down to 160028 peaks (the number of peaks in HCASMC ATACseq)

library(data.table)
library(dplyr)
library(dtplyr)
library(gtools)
library(stringr)

# Constant: 
OUT_DIR='../processed_data/hcasmc_specific_open_chromatin/encode_plus_hcasmc_filt/'



# Functions:
filter_metadata=function(x){
	# filter out ENCODE samples with 'AUDIT ERROR' or 'AUDIT NOT COMPLIANT' flag.
	y=x%>%filter(Assembly=='hg19',`Audit ERROR`=='',`Audit NOT_COMPLIANT`=='')
	y[,max_size:=max(Size),by=`Biosample term name`]
	z=y[Size==max_size,.(`File accession`,`Biosample term name`,`Biosample type`)]
	z[,`File accession`:=paste0(`File accession`,'.bed.gz')]
	setnames(z,c('file','epigenome','type'))
	return(z)
}


sort_bed_by_chrom_pos=function(x){
	# Sort bed file by chromosome and start position
	x=as.data.frame(x)
	x$sortcol=paste(x[,1],x[,2],sep='_')
	ord=mixedorder(x[,'sortcol'])
	y=x[ord,]
	y$sortcol=NULL
	y=as.data.table(y)
}


get_num_peaks=function(){
	# Get number of HCASMC peaks after keeping only autosomes. 
	command=paste('zcat ../processed_data/hcasmc_specific_open_chromatin/encode_plus_hcasmc/',annotation$file[1],sep='/')
	x=fread(command)
	setnames(x,c('chrom','chromStart','chromEnd','name','score','strand','signalValue','pValue','qValue','peak'))
	y=x[chrom%in%paste0('chr',1:22),]
	return(nrow(y))
}


# Rename sample by cell type: 
annotation=fread('../data/encode/dnase_seq/metadata.tsv')%>%filter(Assembly=='hg19')
annotation=filter_metadata(annotation)
annotation=rbind(data.frame(file='2305_ppr.IDR0.1.filt.narrowPeak.gz',epigenome='HCASMC',type='primary cell'),annotation)



# Process DHS data:
# num_peaks=get_num_peaks()
message(sprintf('%d inputs in total',nrow(annotation)))
# for (i in 1:nrow(annotation)){
# 	input=annotation$file[i]
# 	message(input)
# 	command=paste('zcat ../processed_data/hcasmc_specific_open_chromatin/encode_plus_hcasmc/',input,sep='/')
# 	x=fread(command)
# 	setnames(x,c('chrom','chromStart','chromEnd','name','score','strand','signalValue','pValue','qValue','peak'))
# 	y=x[chrom%in%paste0('chr',1:22),] # keeping only autosomes.
# 	y=y%>%arrange(desc(signalValue)) # sort peaks by signalValue
# 	if (nrow(y)>=num_peaks){
# 		z=y[1:num_peaks,] # take the top 50,000 peaks 
# 	} else {
# 		message(sprintf('%s has less than %d peaks',input,num_peaks))
# 		z=y
# 	}
# 	z=sort_bed_by_chrom_pos(z)
# 	output=paste0(OUT_DIR,'/',str_replace_all(annotation$epigenome[i],' ','_'),'.bed')
# 	fwrite(z,output,sep='\t',col.names=F)
# }

for (i in 1:nrow(annotation)){
	input=annotation$file[i]
	output=paste0(OUT_DIR,'/',str_replace_all(annotation$epigenome[i],' ','_'),'.bed')
	message(input)
	command=sprintf('zcat ../processed_data/hcasmc_specific_open_chromatin/encode_plus_hcasmc/%s | sort -k1,1 -k2,2n | bedtools merge -c 7 -o sum -i - | grep -v -e chrX -e chrY -e chrM > %s',input,output)
	system(command)
}