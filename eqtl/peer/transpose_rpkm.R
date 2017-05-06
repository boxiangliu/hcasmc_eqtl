#!/usr/bin/env Rscript
# bosh liu
# durga
# 160528
# transpose the rpkm matrix

# command arguments: 
args=commandArgs(T)
input_file=args[1]
output_file=args[2]


# read input:
input=fread(input_file,header=T)


# separate gene name from rpkm: 
row_names=input$Name
stopifnot(length(row_names)>0)
rpkm=as.matrix(input[,-'Name',with=F])
stopifnot(ncol(rpkm)==ncol(input)-1)
rownames(rpkm)=row_names


# transpose:
output=t(rpkm)


# write output:
write.table(output,file=output_file,quote=F,sep='\t',row.names=T,col.names=T)

