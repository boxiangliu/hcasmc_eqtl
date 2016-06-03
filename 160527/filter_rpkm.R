#!/usr/bin/env Rscript
# bosh liu
# durga
# filter rpkm

args=commandArgs(T)
min_rpkm=args[1]
# min_rpkm=0.1
min_num_indv=args[2]
# min_num_indv=10
input_file=args[3]
# input_file='../processed_data/160527/combined.rpkm'
output_file=args[4]
# output_file='../processed_data/160527/combined.filter.rpkm'

# read input:
input=fread(input_file,header=T)


# keep rows with min_rpkm on at least min_num_indv:
num_row=nrow(input)
keep=rep(T,num_row)
for (i_row in 1:num_row){
	keep[i_row]=sum(input[i_row,-1,with=F]>min_rpkm)>=min_num_indv
}
input_keep=input[keep,]


# write filtered input: 
write.table(input_keep,file=output_file,quote=F,sep='\t',row.names=F,col.names=T)
