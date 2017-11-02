library(data.table)
library(xlsx)
library(stringr)

# Variables:
in_dir='../data/encode/dnase_seq/'
in_dir_2007_2012='../data/encode/dnase_seq_2007_2012/'

out_dir='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks/'
out_dir_filt='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_filt/'
out_dir_released='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_released/'
out_dir_released_all_cell_type='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_released_all_cell_type/'
out_dir_adult='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_adult/' # subsetted to adult samples (no fetal, child, postnatal or newborn)
out_dir_adult_filt='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_adult_filt/'
out_dir_2007_2012='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_2007_2012/'
out_dir_2007_2012_noCancer='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_2007_2012_noCancer/'

for (d in c(out_dir,
			out_dir_filt,
			out_dir_released,
			out_dir_released_all_cell_type,
			out_dir_adult,
			out_dir_adult_filt,
			out_dir_2007_2012,
			out_dir_2007_2012_noCancer)){
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


merge_encode_2007_2012_cell_lines=function(metadata,in_dir,out_dir){
	for (t in unique(metadata[,broad_tissue_category])){
		print(sprintf('INFO - %s',t))
		print(nrow(metadata[broad_tissue_category==t,]))
		out_fn=sprintf('%s/%s.merged.bed',out_dir,t)
		if (file.exists(out_fn)){
			print(sprintf('INFO - %s exist. Skipping...',out_fn))
			next()
		}
		container=list()
		for (fn in metadata[broad_tissue_category==t,BED_File]){
			fn=list.files(in_dir,pattern=fn,full.names=T,recursive=T)
			print(sprintf("INFO - reading %s",fn))
			x=fread(sprintf('zcat %s',fn),col.names=c('chr','start','end','name','score','strand','signalValue','pValue','qValue','peak'))
			if (!all(unlist(x[1:10,list((end-start)==150)]))){
				x[,c('start','end'):=list(start+peak-75,start+peak+75)]
			}
			stopifnot(all(unlist(x[1:10,list((end-start)==150)])))
			container[[fn]]=x[,list(chr,start,end)]
		}
		y=unique(Reduce(rbind,container))
		fwrite(y,sprintf('%s/%s.merged.bed',out_dir,t),sep='\t')
	}
	print('Done!')
}

# HCASMC:
in_fn=list.files('../data/atacseq/fbs/',full.names=T,pattern='naive_overlap.narrowPeak.gz',recursive=T)

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
setorder(dhs,chr,start)


for (d in c(out_dir,
			out_dir_filt,
			out_dir_released,
			out_dir_released_all_cell_type,
			out_dir_adult,
			out_dir_adult_filt,
			out_dir_2007_2012,
			out_dir_2007_2012_noCancer)){
	fwrite(dhs,sprintf('%s/HCASMC.merged.bed',d),sep='\t')
}



# ENCODE (all life stages):
metadata=fread('../data/encode/dnase_seq/metadata.tsv')
metadata=metadata[`Biosample type`%in%c('tissue','primary cell')]
merge_encode(metadata)


# ENCODE (all life stages, remove low quality sample):
metadata=fread('../data/encode/dnase_seq/metadata.tsv')
metadata=metadata[`Biosample type`%in%c('tissue','primary cell')&`Audit ERROR`=='']
merge_encode(metadata)


# ENCODE (all released samples):
metadata=fread('../data/encode/dnase_seq/metadata.tsv')
metadata=metadata[`Biosample type`%in%c('tissue','primary cell')&`File Status`=='released']
# length(unique(metadata$`Biosample term name`)) # 128
# length(metadata$`Biosample term name`) # 593
merge_encode(metadata,in_dir,out_dir_released)


# ENCODE (all released samples for all cell types):
metadata=fread('../data/encode/dnase_seq/metadata.tsv')
metadata=metadata[`File Status`=='released']
metadata[,`Biosample term name`:=str_replace_all(`Biosample term name`,'\\/','_')]
merge_encode(metadata,in_dir,out_dir_released_all_cell_type)


# ENCODE (life stages=adult):
metadata=fread('../data/encode/dnase_seq/metadata.tsv')
metadata=metadata[`Biosample type`%in%c('tissue','primary cell')&`Biosample life stage`=='adult']
merge_encode(metadata)


# ENCODE (life stages=adult, remove low quality sample):
metadata=fread('../data/encode/dnase_seq/metadata.tsv')
metadata=metadata[`Biosample type`%in%c('tissue','primary cell')&`Biosample life stage`=='adult'&`Audit ERROR`=='']
merge_encode(metadata)





# ENCODE 125 Uniformly processed cell lines, categorized into tissue groups:
metadata=read.xlsx('gwas_atacseq_overlap/papers/GREGOR_SuppTable1.xlsx',sheetIndex=1,startRow=3, header=TRUE,colClasses="character")
metadata=as.data.table(apply(metadata,2,as.character))
merge_encode_2007_2012_cell_lines(metadata,in_dir_2007_2012,out_dir_2007_2012)



# ENCODE 125 Uniformly processed cell lines (minus cancer lines), categorized into tissue groups:
encode_human_cell_type_fn='../data/encode/dnase_seq_2007_2012/controlled_vocabulary/human_cell_types.txt'
encode_human_cell_type=fread(encode_human_cell_type_fn)
encode_human_cell_type[term=='Ishikawa',karyotype:='cancer']


sample_annotation_fn=list.files('../data/encode/dnase_seq_2007_2012/',pattern='files.txt',full.names=TRUE,recursive=TRUE)

container=list()
for (f in sample_annotation_fn){
	container[[f]]=fread(f,header=FALSE,col.names=c('filename','annotation'))
}
sample_annotation=Reduce(rbind,container)

sample_annotation=sample_annotation[str_detect(filename,'narrowPeak')]
sample_annotation[,filename:=str_replace(filename,'.gz','')]
sample_annotation[,c('dataType',
					 'cell',
					 'treatment',
					 'type'):=parse_annotation(annotation)]
sample_annotation[,annotation:=NULL]

sample_annotation=merge(sample_annotation,encode_human_cell_type[,list(term,karyotype,tissue)],by.x='cell',by.y='term')
sample_annotation_noCancer=sample_annotation[karyotype!='cancer']

metadata=read.xlsx('gwas_atacseq_overlap/papers/GREGOR_SuppTable1.xlsx',sheetIndex=1,startRow=3, header=TRUE,colClasses="character")
metadata=as.data.table(apply(metadata,2,as.character))

metadata_noCancer=merge(metadata,sample_annotation_noCancer,by.x='BED_File',by.y='filename')

merge_encode_2007_2012_cell_lines(metadata_noCancer,in_dir_2007_2012,out_dir_2007_2012_noCancer)
