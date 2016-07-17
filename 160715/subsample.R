# subsample down to specific input file columns based on column names.

# function:
source('/srv/persistent/bliu2/HCASMC_eQTL/scripts/utils.R')


# load read counts: 
input_file='/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/v6p/v6p_All_Tissues_read_counts_FOR_QC_ONLY.gct'
input=fread(input_file,skip=2,header=T)


# load sample files: 
sample_files=list.files('/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/v6p/subsampling/',pattern='*.52.txt',full.names=T,recursive=T)


# reformat input and decompose it into the count table, 
# the row_data and the column data
res=decompose_gct(input)
count=res$count
col_data=res$col_data
row_data=res$row_data


# loop over sample_files:
for (sample_file in sample_files){
	message(sample_file)


	# read counts: 
	sample=fread(sample_file,header=F)


	# subset samples: 
	sample_count=count[colnames(count)%in%sample$V1]


	# write output:
	sample_count=data.table(row_data,sample_count)
	output_file=str_replace(sample_file,'txt','count')
	write.table(sample_count,file=output_file,quote=F,sep='\t',row.names=F)
}


