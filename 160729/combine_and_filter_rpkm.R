#!/usr/bin/env Rscript 
# bosh liu
# durga
# combine RPKM across all samples listed in "input_list_file", and 
# only keep genes with expression > "min_rpkm" for > "min_prop"
# individuals. 


# command args: 
args=commandArgs(T)
input_list_file=args[1]
out_prefix=args[2]
min_rpkm=as.numeric(args[3])
min_prop=as.numeric(args[4])

# input_list_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/sample_list.txt.head'
# out_prefix='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined'
# min_rpkm=0.1
# min_prop=0.8

# read input list:
input_list=scan(input_list_file,what=character()) 


# read and combine rpkm files:
combined_rpkm=data.frame()
col_data=data.frame()
for (input in input_list){
	message('reading ',input)
	tissue=input %>% basename() %>% str_replace('_Analysis.v6p.FOR_QC_ONLY.rpkm','')

	# read a rpkm file: 
	rpkm=fread(input,header=T)
	stopifnot(length(unique(rpkm$Name))==length(rpkm$Name)) # sanity check

	# if this rpkm file is the first being read: 
	if (ncol(combined_rpkm)==0){
		row_data=rpkm$Name
		rpkm[,Name:=NULL] 
		combined_rpkm=rpkm
		col_data=data.frame(sample=colnames(rpkm),tissue=tissue)
	} 
	# if this rpkm file is not the first being read: 
	else {
		# make the row orders of this rpkm dataframe consistent with previous rpkm dataframes: 
		if (!all(rpkm$Name==row_data)) {
			stopifnot(setequal(rpkm$Name,row_data))
			idx=match(row_data,rpkm$Name)
			rpkm=rpkm[idx,]
		}
		# merge the rpkm dataframe with previous rpkm dataframes:
		stopifnot(rpkm$Name==row_data)
		rpkm[,Name:=NULL]
		combined_rpkm=cbind(combined_rpkm,rpkm)
		col_data=rbind(col_data, data.frame(sample=colnames(rpkm),tissue=tissue))
	}
}

# filter for genes with greater than min. rpkm in min. proportion of samples :
message('filtering...')
pass=combined_rpkm>min_rpkm
n_sample=ncol(pass)
keep=rowSums(pass)>(min_prop*n_sample)
combined_rpkm=combined_rpkm[keep,]
row_data=row_data[keep]


# log2(x+1) transform:
message('log transforming...')
combined_rpkm=log2(combined_rpkm+1)


# add gene id: 
combined_rpkm$Name=row_data
setcolorder(combined_rpkm,c(ncol(combined_rpkm),seq(ncol(combined_rpkm)-1)))


# write output: 
out_rpkm=paste0(out_prefix,'.rpkm')
out_col_data=paste0(out_prefix,'.col')
message('writing output...')
write.table(combined_rpkm,out_rpkm,sep="\t",quote=F,col.names=T,row.names=F)
write.table(col_data,out_col_data,sep="\t",quote=F,col.names=T,row.names=F)
