library(data.table)
library(xlsx)

# Variables:
in_dir='../data/encode/dnase_seq/'
in_dir_2007_2012='../data/encode/dnase_seq_2007_2012/'

out_dir='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks/'
out_dir_adult='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_adult/' # subsetted to adult samples (no fetal, child, postnatal or newborn)
out_dir_adult_filt='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_adult_filt/'
out_dir_2007_2012='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_2007_2012/'


if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
if (!dir.exists(out_dir_adult)) {dir.create(out_dir_adult,recursive=TRUE)}
if (!dir.exists(out_dir_adult_filt)) {dir.create(out_dir_adult_filt,recursive=TRUE)}
if (!dir.exists(out_dir_2007_2012)) {dir.create(out_dir_2007_2012,recursive=TRUE)}

# HCASMC:
in_fn=list.files('../data/atacseq/tmp_from_bliu2_atac_hcasmc/fbs/',full.names=T,pattern='naive_overlap.narrowPeak.gz',recursive=T)

container=list(length(in_fn))
for (i in 1:length(in_fn)){
	dhs=fread(sprintf('zcat %s',in_fn[[i]]),col.names=c('chr','start','end','name','score','strand','signalValue','pValue','qValue','peak'))
	container[[i]]=dhs
}

dhs=Reduce(rbind,container)
dhs=dhs[,list(chr,start,end,peak)]
dhs[,c('start','end'):=list(start+peak-75,start+peak+75)]
dhs[,peak:=NULL]
dhs=unique(dhs)
fwrite(dhs,sprintf('%s/HCASMC.merged.bed',out_dir),sep='\t')
fwrite(dhs,sprintf('%s/HCASMC.merged.bed',out_dir_adult),sep='\t')
fwrite(dhs,sprintf('%s/HCASMC.merged.bed',out_dir_adult_filt),sep='\t')
fwrite(dhs,sprintf('%s/HCASMC.merged.bed',out_dir_2007_2012),sep='\t')


# ENCODE (all life stages):
metadata=fread('../data/encode/dnase_seq/metadata.tsv')
metadata=metadata[`Biosample type`%in%c('tissue','primary cell')]
biosample=unique(metadata$`Biosample term name`)
for (b in biosample){
	print(sprintf('INFO - %s',b))
	tmp=metadata[`Biosample term name`==b,list(`File accession`,`Biosample term name`,peak_type)]

	container=list()
	for (i in 1:nrow(tmp)){
		x=fread(sprintf('zcat %s/%s.bed.gz',in_dir,tmp$`File accession`[i]),col.names=c('chr','start','end','name','score','strand','signalValue','pValue','qValue','peak'))

		if (tmp$peak_type[i]=='raw'){
			x[,c('start','end'):=list(start+peak-75,start+peak+75)]
		}
		container[[i]]=x[,list(chr,start,end)]
	}
	y=unique(Reduce(rbind,container))
	fwrite(y,sprintf('%s/%s.merged.bed',out_dir,b),sep='\t')
}



# ENCODE (life stages=adult):
metadata=fread('../data/encode/dnase_seq/metadata.tsv')
metadata=metadata[`Biosample type`%in%c('tissue','primary cell')&`Biosample life stage`=='adult']
biosample=unique(metadata$`Biosample term name`)

for (b in biosample){
	print(sprintf('INFO - %s',b))
	tmp=metadata[`Biosample term name`==b,list(`File accession`,`Biosample term name`,peak_type)]

	container=list()
	for (i in 1:nrow(tmp)){
		x=fread(sprintf('zcat %s/%s.bed.gz',in_dir,tmp$`File accession`[i]),col.names=c('chr','start','end','name','score','strand','signalValue','pValue','qValue','peak'))

		if (tmp$peak_type[i]=='raw'){
			x[,c('start','end'):=list(start+peak-75,start+peak+75)]
		}
		container[[i]]=x[,list(chr,start,end)]
	}
	y=unique(Reduce(rbind,container))
	fwrite(y,sprintf('%s/%s.merged.bed',out_dir_adult,b),sep='\t')
}


# ENCODE (life stages=adult, remove low quality sample):
metadata=fread('../data/encode/dnase_seq/metadata.tsv')
metadata=metadata[`Biosample type`%in%c('tissue','primary cell')&`Biosample life stage`=='adult'&`Audit ERROR`=='']
biosample=unique(metadata$`Biosample term name`)

for (b in biosample){
	print(sprintf('INFO - %s',b))
	tmp=metadata[`Biosample term name`==b,list(`File accession`,`Biosample term name`,peak_type)]

	container=list()
	for (i in 1:nrow(tmp)){
		x=fread(sprintf('zcat %s/%s.bed.gz',in_dir,tmp$`File accession`[i]),col.names=c('chr','start','end','name','score','strand','signalValue','pValue','qValue','peak'))

		if (tmp$peak_type[i]=='raw'){
			x[,c('start','end'):=list(start+peak-75,start+peak+75)]
		}
		container[[i]]=x[,list(chr,start,end)]
	}
	y=unique(Reduce(rbind,container))
	fwrite(y,sprintf('%s/%s.merged.bed',out_dir_adult_filt,b),sep='\t')
}

# ENCODE 125 Uniformly processed cell lines, categorized into tissue groups:
metadata=read.xlsx('gwas_atacseq_overlap/papers/GREGOR_SuppTable1.xlsx',sheetIndex=1,startRow=3, header=TRUE,colClasses="character")
metadata=as.data.table(apply(metadata,2,as.character))
for (tissue in unique(metadata[,broad_tissue_category])){
	print(sprintf('INFO - %s',tissue))
	out_fn=sprintf('%s/%s.merged.bed',out_dir_2007_2012,tissue)
	if (file.exists(out_fn)){
		print(sprintf('INFO - %s exist. Skipping...',out_fn))
		next()
	}
	container=list()
	for (fn in metadata[broad_tissue_category==tissue,BED_File]){
		fn=list.files(in_dir_2007_2012,pattern=fn,full.names=T,recursive=T)
		print(sprintf("INFO - reading %s",fn))
		x=fread(sprintf('zcat %s',fn),col.names=c('chr','start','end','name','score','strand','signalValue','pValue','qValue','peak'))
		if (!all(unlist(x[1:10,list((end-start)==150)]))){
			x[,c('start','end'):=list(start+peak-75,start+peak+75)]
		}
		stopifnot(all(unlist(x[1:10,list((end-start)==150)])))
		container[[fn]]=x[,list(chr,start,end)]
	}
	y=unique(Reduce(rbind,container))
	fwrite(y,sprintf('%s/%s.merged.bed',out_dir_2007_2012,tissue),sep='\t')
}