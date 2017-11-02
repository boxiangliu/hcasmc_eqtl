library(data.table)
library(XLConnect)
library(NOISeq)
library(foreach)
library(doMC)
registerDoMC(6)

# Variables:
hcasmc_rpkm_fn='/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/rpkm/combined.rpkm'
hcasmc_covariate_fn='/srv/persistent/bliu2/HCASMC_eQTL/data/sample_info/sample_info.xlsx'

subject_dir='../processed_data/differential_expression/sample_individuals/'
gtex_covariate_fn='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/sample_annotations/GTEx_Analysis_2015-01-12_Annotations_SubjectPhenotypesDS.txt'
gtex_rpkm_dir='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6p_rpkm/'

out_dir='../processed_data/differential_expression/noiseq/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

tissue_to_id_map=fread('/srv/persistent/bliu2/HCASMC_eQTL/data//gtex/gtex_tissue_colors.txt',select=c(1,3))
coding_and_lncRNA=get_coding_and_lncRNA()
tissue_list=c("Artery - Aorta","Artery - Coronary","Artery - Tibial", 'Heart - Atrial Appendage', 'Heart - Left Ventricle','Cells - Transformed fibroblasts')



# Functions:
get_coding_and_lncRNA=function(){
	gencode=fread('/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/gencode.v19.genes.v6p.hg19.bed',select=c(5,6,7),col.names=c('gene_id','gene_name','type'))
	coding_and_lncRNA=gencode[type%in%c('lincRNA','protein_coding'),gene_id]
	return(coding_and_lncRNA)
}

read_hcasmc_rpkm=function(fn,subset){
	x=fread(fn,head=TRUE)
	x=x[Name%in%subset]
	x_mat=as.matrix(x[,2:ncol(x)])
	rownames(x_mat)=x$Name
	return(x_mat)
}


read_hcasmc_covariate=function(hcasmc_covariate_fn){
	hcasmc_covariate=readWorksheet(loadWorkbook(hcasmc_covariate_fn),sheet=5)
	setDT(hcasmc_covariate)
	hcasmc_covariate=hcasmc_covariate[,list(SUBJID=DNA,sex=Sex,age=Age_Imputed,race=Genomic_Ethnicity)]
	return(hcasmc_covariate)
}


merge_count_table=function(hcasmc_read_count_mat,gtex_read_count_downsample_mat){
	row_idx=match(rownames(hcasmc_read_count_mat),rownames(gtex_read_count_downsample_mat))
	gtex_read_count_downsample_mat_reordered=gtex_read_count_downsample_mat[row_idx,]
	stopifnot(rownames(hcasmc_read_count_mat)==rownames(gtex_read_count_downsample_mat_reordered))
	read_count_mat=cbind(hcasmc_read_count_mat,gtex_read_count_downsample_mat_reordered)
	return(read_count_mat)
}


get_subject_list=function(suject_dir,tissue){
	subject_fn=sprintf('%s/%s.txt',subject_dir,tissue)
	subject_list=fread(subject_fn,header=FALSE,col.names='subject')
	return(subject_list)
}


read_gtex_rpkm=function(dir,tissue_id,subject_list){
	fn=sprintf('%s/%s.rpkm.txt',dir,tissue_id)
	x=fread(fn)
	x=x[Gene%in%coding_and_lncRNA]
	x_mat=as.matrix(x[,unlist(subject_list),with=F])
	rownames(x_mat)=x[,Gene]
	return(x_mat)
}


read_gtex_covariate=function(gtex_covariate_fn,subject_list){
	gtex_covariate=fread(gtex_covariate_fn)[,list(SUBJID,age=AGE,sex=GENDER,race=RACE,ethnicity=ETHNCTY)]
	race_numeric_to_letter=c('Asian','AA','Caucasian','Indian','Unknown')
	gtex_covariate[,race:=race_numeric_to_letter[race]]
	gtex_covariate[,race:=ifelse(ethnicity==1,'Hispanic',race)]
	gtex_covariate[,ethnicity:=NULL]
	sex_numeric_to_letter=c('M','F')
	gtex_covariate[,sex:=sex_numeric_to_letter[sex]]
	tissue_covariate=gtex_covariate[SUBJID%in%unlist(subject_list)]
	return(tissue_covariate)
}


merge_covariates=function(hcasmc_covariate,tissue_covariate,tissue1='HCASMC',tissue2=tissue){
	hcasmc_covariate$tissue=tissue1
	tissue_covariate$tissue=tissue2
	covariate=rbind(hcasmc_covariate,tissue_covariate)
	return(covariate)
}


reorder_covariate_by_read_count_matrix=function(covariate,read_count){
	covariate[match(colnames(read_count),covariate$SUBJID)]
}



construct_biological_annotation=function(subset){
	gencode=fread('/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/gencode.v19.genes.v6p.hg19.bed',col.names=c('Chr','GeneStart','GeneEnd','Strand','gene_id','gene_name','type'))
	gencode=gencode[match(subset,gencode$gene_id)]


	chroms=as.data.frame(gencode[,list(Chr,GeneStart,GeneEnd)])
	rownames(chroms)=gencode$gene_id
	biotypes=gencode$type
	names(biotypes)=gencode$gene_id
	length=gencode[,abs(GeneEnd-GeneStart)]
	names(length)=gencode$gene_id

	temp=fread('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/rasqual/expression/gcc.exon.txt',col.names=c('gene_id','gc'))
	temp=temp[match(subset,temp$gene_id)]

	gc=temp$gc
	names(gc)=temp$gene_id
	return(list(length=length,gc=gc,biotypes=biotypes,chroms=chroms))
}


get_noiseqbio_result=function(x){
	results=as.data.table(nsb@results,keep.rowname=TRUE)
}



# Main:
# Read HCASMC RPKM and covariate:
hcasmc_rpkm=read_hcasmc_rpkm(hcasmc_rpkm_fn,subset=coding_and_lncRNA)
hcasmc_covariate=read_hcasmc_covariate(hcasmc_covariate_fn)


container=foreach(tissue=tissue_list,.final=function(x) setNames(x,tissue_list))%dopar%{
	message(tissue)

	# Read GTEx RPKM and covariate:
	tissue_id=tissue_to_id_map[tissue_site_detail==tissue,tissue_site_detail_id]
	subject_list=get_subject_list(subject_dir,tissue)
	gtex_rpkm=read_gtex_rpkm(gtex_rpkm_dir,tissue_id,subject_list)
	gtex_covariate=read_gtex_covariate(gtex_covariate_fn,subject_list)


	# Merge RPKM and covariate:
	rpkm=merge_count_table(hcasmc_rpkm,gtex_rpkm)
	covariate=merge_covariates(hcasmc_covariate,gtex_covariate,tissue1='HCASMC',tissue2=tissue_id)
	covariate=reorder_covariate_by_read_count_matrix(covariate,rpkm)


	# Run NOISeq:
	bio_anno=construct_biological_annotation(rownames(rpkm))
	data=readData(data=rpkm,length=bio_anno$length,
		gc=bio_anno$gc,biotype=bio_anno$biotypes,
		chromosome=bio_anno$chroms,factors=as.data.frame(covariate[,list(sex,race,tissue)]))
	nsb = noiseqbio(data, k = 0.5, norm = "n", 
		factor = "tissue", r = 20,adj = 1.5, plot = FALSE, a0per = 0.9, 
		random.seed = 12345, filter = 2)

	return(nsb)
}

# Save result:
saveRDS(container,sprintf('%s/noiseq.rds',out_dir))