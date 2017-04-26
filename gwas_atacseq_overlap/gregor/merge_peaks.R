library(data.table)

# Variables:
out_dir='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks/'
out_dir_adult='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_adult/' # subsetted to adult samples (no fetal, child, postnatal or newborn)
in_dir='../data/encode/dnase_seq/'

if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
if (!dir.exists(out_dir_adult)) {dir.create(out_dir_adult,recursive=TRUE)}


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
