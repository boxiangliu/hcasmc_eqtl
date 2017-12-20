#!/usr/bin/env Rscript 
# boxiang liu
# durga
# combine RPKM across 53 GTEx tissues and HCASMC
library(data.table)

# read and combine rpkm files:
gtex_rpkm_fn='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6p_All_Tissues_gene_rpkm_FOR_QC_ONLY.gct'
gtex_sample_annotation_fn='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/sample_annotations/GTEx_Analysis_2015-01-12_Annotations_SampleAttributesDS.txt'
hcasmc_rpkm_fn='../processed_data/rnaseq/preprocess/combine_rpkm/combined.rpkm'
out_dir='../processed_data/mds/combine_hcasmc_and_gtex_rpkm/'
rpkm_threshold=0.1
num_indv_threshold=10
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}


# Function:
get_col_data=function(combined_rpkm,gtex_sample_annotation){
	col_data=data.table(sample_id=colnames(combined_rpkm)[2:ncol(combined_rpkm)])
	tmp=data.table(sample_id=colnames(hcasmc)[2:ncol(hcasmc)],tissue='HCASMC')
	sample_annotation=rbind(gtex_sample_annotation,tmp)
	col_data=merge(col_data,sample_annotation,by='sample_id')
	return(col_data)
}


# Read GTEx RPKMs:
gtex=fread(gtex_rpkm_fn)
gtex_sample_annotation=fread(gtex_sample_annotation_fn)[,list(sample_id=SAMPID,tissue=SMTSD)]

# Read HCASMC RPKM:
hcasmc=fread(hcasmc_rpkm_fn,header=TRUE)

# Combined GTEx and HCASMC RPKMs: 
combined_rpkm=merge(gtex,hcasmc,by='Name')

# filter for genes with rpkm>0.1 in > 10 individuals:
message('filtering...')
pass=combined_rpkm[,3:ncol(combined_rpkm)]>rpkm_threshold
keep=rowSums(pass)>num_indv_threshold
combined_rpkm=combined_rpkm[keep,]
row_data=combined_rpkm[,list(Name,Description)]
fwrite(row_data,sprintf('%s/combined.row',out_dir),sep='\t')

# log2(x+2) transform:
message('log transforming...')
combined_rpkm=log2(combined_rpkm[,3:ncol(combined_rpkm)]+2)

# add gene id: 
combined_rpkm$Name=row_data$Name
setcolorder(combined_rpkm,c(ncol(combined_rpkm),seq(ncol(combined_rpkm)-1)))

# write output:
message('writing output...')
fwrite(combined_rpkm,sprintf('%s/combined.rpkm',out_dir),sep='\t')
col_data=get_col_data(combined_rpkm,gtex_sample_annotation)
fwrite(col_data,sprintf('%s/combined.col',out_dir),sep='\t')
