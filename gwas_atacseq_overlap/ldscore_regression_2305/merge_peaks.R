library(data.table)
library(xlsx)
library(stringr)

# Variables:
in_dir='../data/encode/dnase_seq/'

out_dir_released='../processed_data/gwas_atacseq_overlap/ldscore_regression_2305/merge_peaks_released/'
out_dir_released_all_cell_type='../processed_data/gwas_atacseq_overlap/ldscore_regression_2305/merge_peaks_released_all_cell_type/'

for (d in c(out_dir_released,out_dir_released_all_cell_type)){
	if (!dir.exists(d)) {dir.create(d,recursive=TRUE)}
}

# Functions:
parse_annotation=function(annotation){
	dataType=str_extract(annotation,'(?<=dataType=)(.+?)(?=;)')
	cell=str_extract(annotation,'(?<=cell=)(.+?)(?=;)')
	treatment=str_extract(annotation,'(?<=treatment=)(.+?)(?=;)')
	type=str_extract(annotation,'(?<=type=)(.+?)(?=;)')
	list(dataType,
		 cell,
		 treatment,
		 type)
}


merge_encode=function(metadata,in_dir,out_dir){
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
}

# HCASMC:
in_fn='../data/atacseq/fbs/2305/out/peak/macs2/overlap/CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.nodup.tn5_pooled.pf.pval0.1.500K.naive_overlap.narrowPeak.gz'
dhs=fread(sprintf('zcat %s',in_fn),col.names=c('chr','start','end','name','score','strand','signalValue','pValue','qValue','peak'))
dhs=dhs[,list(chr,start,end,peak)]
dhs[,c('start','end'):=list(start+peak-75,start+peak+75)]
dhs[,peak:=NULL]
dhs=unique(dhs)
setorder(dhs,chr,start)


for (d in c(out_dir_released,out_dir_released_all_cell_type)){
	fwrite(dhs,sprintf('%s/HCASMC.merged.bed',d),sep='\t')
}


# ENCODE (all released samples):
metadata=fread('../data/encode/dnase_seq/metadata.tsv')
metadata=metadata[`Biosample type`%in%c('tissue','primary cell')&`File Status`=='released']
merge_encode(metadata,in_dir,out_dir_released)


# ENCODE (all released samples for all cell types):
metadata=fread('../data/encode/dnase_seq/metadata.tsv')
metadata=metadata[`File Status`=='released']
metadata[,`Biosample term name`:=str_replace_all(`Biosample term name`,'\\/','_')]
merge_encode(metadata,in_dir,out_dir_released_all_cell_type)