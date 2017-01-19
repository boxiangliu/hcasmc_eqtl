library(data.table)
library(stringr)

# functions: 
split_fid=function(fid,sep="_"){
	x=str_split_fixed(fid,pattern=sep,2)
	colnames(x)=c('gene_id','gene_name')
	return(x)
}

# variables:
in_file='../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_variants.allFeatures.bed'
out_file='../processed_data/mpra/forNathan/eCAVIAR_variants.bed'


# main:
variants=fread(in_file)
variants=cbind(variants,split_fid(variants$fid))
columns=c('chr','start','end','ref','alt','gene_id','rsid','clpp','gwasSet','eqtlSet')
output=variants[,columns,with=F]
if(!dir.exists(dirname(out_file))){dir.create(dirname(out_file))}
write.table(output,out_file,sep='\t',col.names=T,row.names=F,quote=F)



