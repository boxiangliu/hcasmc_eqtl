# Select ENCODE cell lines with no "AUDIT NOT COMPLIANT" or "AUDIT ERROR" flag.
# Subsample each selected sample down to 160028 peaks (the number of peaks in HCASMC ATACseq)

library(data.table)
library(dplyr)
library(dtplyr)
library(gtools)
library(stringr)

# Constant: 
IN_DIR='../data/encode/dnase_seq/'
OUT_DIR='../processed_data/hcasmc_specific_open_chromatin2/encode_plus_hcasmc_filt/'
FIG_DIR='../figures/hcasmc_specific_open_chromatin2/process_encode_data/'
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR)

# Functions:
add_metadata=function(metadata){
	if ('peak_type' %in% colnames(metadata) & 'n_peaks' %in% colnames(metadata) & 'window' %in% colnames(metadata)){
		return(metadata)
	}
	fn_ls=metadata[,`File accession`]
	n_peaks=integer(length(fn_ls))
	peak_type=integer(length(fn_ls))
	window_size=integer(length(fn_ls))
	for (i in seq_along(fn_ls)){
		fn=fn_ls[i]
		print(sprintf('%s: %s',i,fn))
		x=fread(sprintf('zcat %s/%s.bed.gz',IN_DIR,fn))
		n_peaks[i]=nrow(x)
		peak_type[i]=ifelse(all(unlist(x[1:20,10])==unlist(x[1,10])),'summit_extend','raw')
		window_size[i]=ifelse(peak_type[i]=='summit_extend',unlist(x[1,3])-unlist(x[1,2]),NaN)
	}
	metadata$n_peaks=n_peaks
	metadata$peak_type=peak_type
	metadata$window=window_size
	return(metadata)
}


filter_metadata=function(x){
	y=x%>%filter(Assembly=='hg19',`Biosample type`%in%c('tissue','primary cell'),`Biosample life stage`%in%c('adult','fetal'),`File Status`=='released')
	setDT(y)
	y[,max_peak:=max(n_peaks),by="Biosample term name"]
	z=y[max_peak==n_peaks,.(`File accession`,`Biosample term name`,`Biosample type`,`Biosample life stage`,peak_type)]
	z[,`File accession`:=paste0(`File accession`,'.bed.gz')]
	setnames(z,c('file','epigenome','type','life_stage','peak_type'))
	return(z)
}


# Rename sample by cell type: 
annotation=fread('../data/encode/dnase_seq/metadata.tsv')%>%filter(Assembly=='hg19')
annotation=add_metadata(annotation)


# Output modified metadata:
write.table(annotation,'../data/encode/dnase_seq/metadata.tsv',quote=F,sep='\t',col.names=T,row.names=F)


# Filter metadata:
annotation=filter_metadata(annotation)
annotation[,file:=sprintf('%s/%s',IN_DIR,file)]
annotation=rbind(data.frame(file='../data/atacseq/fbs/2305/out/peak/idr/optimal_set/2305_ppr.IDR0.1.filt.narrowPeak.gz',epigenome='HCASMC',type='primary cell',life_stage='adult',peak_type='raw'),annotation)



# Process DHS data:
message(sprintf('%d inputs in total',nrow(annotation)))
for (i in 1:nrow(annotation)){
	peak_type=annotation$peak_type[i]
	input=annotation$file[i]
	output=paste0(OUT_DIR,'/',str_replace_all(annotation$epigenome[i],' ','_'),'.bed')
	message(input)
	is_raw=FALSE
	if (peak_type=='raw'){
		is_raw=TRUE
		x=fread(sprintf("zcat %s",input))
		setnames(x,c('chrom','chromStart','chromEnd','name','score','strand','signalValue','pValue','qValue','peak'))
		x[,summit:=chromStart+peak]
		x[,chromStart:=summit-75]
		x[,chromEnd:=summit+75]
		input=tempfile()
		fwrite(x,input,col.names=F,row.names=F,sep='\t')
		system(sprintf('bgzip %s',input))
		input=paste0(input,'.gz')
	}
	command=sprintf('zcat %s | grep -v -e chrX -e chrY -e chrM | sort -k1,1 -k2,2n | bedtools merge -c 7 -o sum -i - > %s',input,output)
	system(command)
	if (is_raw) file.remove(input)
}