#!/usr/bin/env Rscript
# bosh liu
# durga
# quantile normalize rpkm for each sample, and inverse rank normalization for each genes


# library:
library(preprocessCore)
source('utils.R')

# command args: 
args=commandArgs(T)
input_file=args[1]
output_file=args[2]

# read input:
input=fread(input_file,header=T)


# separate input into row name and rpkm values:
row_names=input$Name
rpkm=as.matrix(input[,-"Name",with=F])


# quantile normalize each sample: 
rpkm_norm=normalize.quantiles(rpkm,copy=T)


# inverse rank normalize each gene:
num_row=nrow(rpkm_norm)
rpkm_norm2=rpkm_norm
for (i_row in seq(num_row)){
	x=rpkm_norm[i_row,]
	x_norm=getInverseNormal(x)
	rpkm_norm2[i_row,]=x_norm
}


# write output:
colnames(rpkm_norm2)=colnames(rpkm)
output=data.table(Name=row_names, rpkm_norm2)
write.table(output, file=output_file, quote=F, sep='\t', row.names=F, col.names=T)
