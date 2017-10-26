library(data.table)
library(stringr)

# Variables:
sample_attributes_fn='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/sample_annotations/GTEx_Analysis_2015-01-12_Annotations_SampleAttributesDS.txt'
out_dir='../processed_data/differential_expression/sample_individuals/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}


# Read sample attributes: 
sample_attributes=fread(sample_attributes_fn)[SMAFRZE=='USE ME',][str_detect(SAMPID,'GTEX'),]
stopifnot(all(sample_attributes[,SMAFRZE=='USE ME']))
sample_attributes[,SUBJID:=paste0('GTEX-',str_split_fixed(SAMPID,'-',5)[,2])]


# Calculate sample size per tissue:
tissue_list=c("Artery - Aorta","Artery - Coronary","Artery - Tibial", 'Heart - Atrial Appendage', 'Heart - Left Ventricle','Cells - Transformed fibroblasts')

tissue_df=sample_attributes[SMTSD%in%tissue_list]
tissue_sample_size=tissue_df[,list(size=length(SUBJID)),by='SMTSD']


# Rank tissue by # of individuals: 
setorder(tissue_sample_size,size)


# Sample individuals: 
downsample_size=52
subject_id_list=list()
set.seed(42)

for (tissue in tissue_sample_size[,SMTSD]){

	subjects=tissue_df[SMTSD==tissue,SUBJID]
	unselected_subjects=subjects[!subjects%in%unlist(subject_id_list)]
	subject_id_list[[tissue]]=sample(unselected_subjects,downsample_size)

	if (length(subject_id_list[[tissue]])<downsample_size){
		message('insufficient sample size')
		selected_subjects=subjects[subjects%in%unlist(subject_id_list)]
		temp=sample(selected_subjects,downsample_size-subject_id_list[[tissue]])
		subject_id_list[[tissue]]=c(subject_id_list[[tissue]],temp)
	}
}


# Output sampled individuals:
for (tissue in names(subject_id_list)){
	subject_id=subject_id_list[[tissue]]
	writeLines(subject_id,sprintf('%s/%s.txt',out_dir,tissue))
}

