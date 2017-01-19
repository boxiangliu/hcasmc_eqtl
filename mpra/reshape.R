library(reshape2)
in_file=commandArgs(T)[1]
out_file=commandArgs(T)[2]
# in_file='/srv/scratch/bliu2/concat_eCAVIAR_variants//ENSG00000214435.3_AS3M.long'
input=read.table(in_file)
output=dcast(input,V1+V2+V3+V4~V5,value.var='V5')
output=output[,c('V1','V2','V3','V4','gwas','eqtl')]
output$gwas=as.integer(!is.na(output$gwas))
output$eqtl=as.integer(!is.na(output$eqtl))
write.table(output,out_file,quote=F,sep='\t',col.names=F,row.names=F)
