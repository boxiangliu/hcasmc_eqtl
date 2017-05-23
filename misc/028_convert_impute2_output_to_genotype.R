#!/usr/bin/env Rscript
# bosh liu
# 2016/05/05
# durga
# convert impute2 output to genotype (0,1,2)

library(data.table)
source('utils.R')

# read in the first few lines of pre-imputation vcf file: 
vcf_filename="../data/joint/recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.vcf.gzary"
vcf_file=scan(vcf_filename,nmax=30,what='character',sep='\n')


# parse the sample names from the first lines of vcf files:
sample_names=vcf_file[29]
sample_names=str_split(sample_names,pattern='\t')[[1]]
sample_names=sample_names[!(sample_names%in%c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"))]
sample_names=str_split_fixed(sample_names,"_",2)[,1]


# loop over all impute2 output: 
imputation_dir='../data/joint/imputation'
input_filenames=list.files(imputation_dir,pattern='chr.+\\.impute2',full.name=T)
# input_filenames=input_filenames[20:22]


# output_dir
output_dir='../processed_data/028_imputed_genotype/'

for (input_filename in input_filenames){
	# extract chromsome from filename: 
	input_basename=basename(input_filename)
	chr=str_split_fixed(input_basename,"chr|\\.impute2",3)[1,2]


	# create name for output file:
	output=paste0(output_dir,'/chr',chr,'.genotype')
	

	# check if output already exist:
	if (file.exists(output)) {message(output,' already exist. skip.'); next()}


	# print input filename:
	message('input:',input_filename)


	# read impute2 output: 
	impute2=fread(input_filename,header=F)


	# extract genotype:
	genotypes=convImpute2ToGeno(impute2)


	# make SNP ID:
	rownames(genotypes)=impute2[,paste(chr,V3,V4,V5,sep="_")]


	# convert genotypes to data.table:
	genotypes=data.table(genotypes,keep.rownames=T)
	setnames(genotypes,c('SNP',sample_names))


	# save genotypes:
	write.table(genotypes,output,col.names=T,row.names=F,quote=F,sep='\t')
}


