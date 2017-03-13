library(data.table)

# Variables: 
in_dir='../processed_data/hcasmc_specific_open_chromatin/intersect/samplewise_count/'
out_file='../processed_data/hcasmc_specific_open_chromatin/intersect/intersect.count'


# Get all input file names: 
fn_list=list.files(in_dir,pattern='*.count',full.names=T)
x=fread('../processed_data/hcasmc_specific_open_chromatin/intersect/samplewise_count//intersect.IMR-90.count')


# Read and merge input files: 
y=data.table()
for (fn in fn_list){
	print(fn)
	x=fread(fn)
	if (nrow(y)==0){
		y=x
	} else {
		y=merge(y,x,by='id',all=TRUE)
	}
}


# Replace missing values with zeros:
for (j in seq_along(y)){
	set(y,i=which(is.na(y[[j]])),j=j,value=0)
}


# Output:
fwrite(y,out_file,sep='\t',col.names=T,row.names=F)