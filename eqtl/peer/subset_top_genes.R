#!/usr/bin/env Rscript
# bosh liu
# durga
# subset top n genes according to rpkm sum

# command args: 
args=commandArgs(T)
input_file=args[1]
output_file=args[2]
n=args[3]


# read input: 
input=fread(input_file,header=T)


# calculate sum of expression:
rpkm_sum=rowSums(input[,-1,with=F])
rpkm_rank=rank(rpkm_sum)


# subset for top n genes: 
keep=(rpkm_rank<=n)
output=input[keep,]


# write output: 
write.table(output,output_file,quote=F,sep='\t',col.names=T,row.names=F)
