# library: 
library(data.table)
library(dplyr)
library(stringr)
library(cowplot)
library(ggrepel)
library(foreach)

# command args: 
in_file='../processed_data/eqtl/metasoft/metasoft_output/metasoft_output.mcmc.txt'
metasoft_input_fn='../processed_data/eqtl/metasoft/metasoft_input/metasoft_input.txt'
fig_dir='../figures/eqtl/metasoft/tissue_specific_eqtl/'
out_dir='../processed_data/eqtl/metasoft/tissue_specific_eqtl/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

# constants: 
N_TISSUE=45

get_study_name = function(){
	# read tissue name (study name):
	study_name=c("Adipose_Subcutaneous",
	"Adipose_Visceral_Omentum",
	"Adrenal_Gland",
	"Artery_Aorta",
	"Artery_Coronary",
	"Artery_Tibial",
	"Brain_Anterior_cingulate_cortex_BA24",
	"Brain_Caudate_basal_ganglia",
	"Brain_Cerebellar_Hemisphere",
	"Brain_Cerebellum",
	"Brain_Cortex",
	"Brain_Frontal_Cortex_BA9",
	"Brain_Hippocampus",
	"Brain_Hypothalamus",
	"Brain_Nucleus_accumbens_basal_ganglia",
	"Brain_Putamen_basal_ganglia",
	"Breast_Mammary_Tissue",
	"Cells_EBV-transformed_lymphocytes",
	"Cells_Transformed_fibroblasts",
	"Colon_Sigmoid",
	"Colon_Transverse",
	"Esophagus_Gastroesophageal_Junction",
	"Esophagus_Mucosa",
	"Esophagus_Muscularis",
	"HCASMC",
	"Heart_Atrial_Appendage",
	"Heart_Left_Ventricle",
	"Liver",
	"Lung",
	"Muscle_Skeletal",
	"Nerve_Tibial",
	"Ovary",
	"Pancreas",
	"Pituitary",
	"Prostate",
	"Skin_Not_Sun_Exposed_Suprapubic",
	"Skin_Sun_Exposed_Lower_leg",
	"Small_Intestine_Terminal_Ileum",
	"Spleen",
	"Stomach",
	"Testis",
	"Thyroid",
	"Uterus",
	"Vagina",
	"Whole_Blood")
	return(study_name)
}

read_metasoft = function(in_file,study_name = get_study_name()){
	# read metasoft result: 
	header=scan(in_file,what='character',nlines=1,sep='\t')
	metasoft=fread(in_file,skip=1)


	# remove the extra column (issue due to white space)
	metasoft[,V107:=NULL]

	# set metasoft column names:
	col_names=c(paste('pvalue',study_name,sep='_'),paste('mvalue',study_name,sep='_'))
	col_names=c(header[1:16],col_names)
	stopifnot(ncol(metasoft)==length(col_names))
	setnames(metasoft,col_names)

	return(metasoft)
}


read_metasoft_input = function(metasoft_input_fn, study_name = get_study_name()){
	metasoft_input = fread(metasoft_input_fn)
	setDF(metasoft_input)
	rownames(metasoft_input) = unlist(metasoft_input[,1])
	metasoft_input[,1] = NULL
	colnames = foreach (i = study_name,.combine='c')%do%{
		c(paste0(i,'_beta'),paste0(i,'_se'))
	}
	colnames(metasoft_input) = colnames
	return(metasoft_input)
}

subset_mvalue = function(metasoft){
	# subset to mvalues: 
	mvalue=metasoft%>%dplyr::select(contains('mvalue'))%>%as.data.frame()
	stopifnot(ncol(mvalue)==N_TISSUE)
	rownames(mvalue)=metasoft$RSID
	colnames(mvalue)=str_replace(colnames(mvalue),'mvalue_','')
	return(mvalue)
}

subset_pvalue = function(metasoft){
	# subset to pvalues: 
	pvalue=metasoft%>%dplyr::select(contains('pvalue'))%>%dplyr::select(-c(1:5))%>%as.data.frame()
	rownames(pvalue)=metasoft$RSID
	colnames(pvalue)=str_replace(colnames(mvalue),'pvalue_','')
	return(pvalue)
}

subset_beta = function(metasoft_input){
	beta=metasoft_input%>%dplyr::select(contains('beta'))%>%as.data.frame()
	stopifnot(ncol(beta)==N_TISSUE)
	colnames(beta)=str_replace(colnames(beta),'_beta','')
	return(beta)
}


subset_se = function(metasoft_input){
	se=metasoft_input%>%dplyr::select(contains('_se'))%>%as.data.frame()
	stopifnot(ncol(se)==N_TISSUE)
	colnames(se)=str_replace(colnames(se),'_se','')
	return(se)
}

select_tseQTL = function(mvalue,tissue,cutoff = 0.9, cutoff2 = 0.1){
	setDF(mvalue)
	is_eqtl = mvalue[,tissue]>=cutoff

	mvalue[,tissue] = NULL
	indicator = mvalue >= cutoff2

	unique_eqtl = rowSums(indicator,na.rm=TRUE) == 0
	tseQTL_idx = is_eqtl & unique_eqtl
	tseQTL=names(tseQTL_idx)[which(tseQTL_idx)]

	return(tseQTL)
}

subset_protein_coding_and_lncRNA = function(tseQTL){
	annotation = fread(
		'/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.bed',
		select = 5:7,
		col.names = c('gene_id','gene_name','type'))
	split_tseQTL = str_split_fixed(tseQTL,'_',2)
	tseQTL = data.table(
		id = tseQTL,
		gene_id = split_tseQTL[,1],
		snp = split_tseQTL[,2]
		)

	tseQTL = merge(tseQTL, annotation, by = 'gene_id')
	tseQTL = tseQTL[type %in% c('protein_coding','lincRNA')]
	
	return(tseQTL)
}

get_color_map = function(tissue_color_fn='../data/gtex/gtex_tissue_colors.with_hcasmc.txt'){
	tissue_color = fread(tissue_color_fn)
	color_map = tissue_color$tissue_color_hex
	names(color_map) = tissue_color$tissue_site_detail_id
	names(color_map)[names(color_map)=='hcasmc.rpkm'] = 'HCASMC'
	return(color_map)
}

make_forestPM_plot = function(mvalue,pvalue,beta,se,id,color_map,title=''){
	pm_plot_data = data.table(
		tissue = colnames(mvalue),
		mvalue = unlist(mvalue[rownames(mvalue)==id,]),
		pvalue = -log10(unlist(pvalue[rownames(pvalue)==id,]))
		)
	pm_plot_data$label = ifelse(pm_plot_data$tissue=='HCASMC','HCASMC','')
	pm_plot_data = pm_plot_data[!is.na(mvalue)&!is.na(pvalue),]

	p1 = ggplot(pm_plot_data,aes(x = mvalue, y = pvalue, label = label, color = tissue))+
		geom_point(size = 3)+
		xlab('M-value')+
		geom_text(color = 'black', vjust = 2.0, hjust = 1) + 
		ylab('-log10(P-value)')+
		geom_vline(xintercept = 0.9,color = 'red',linetype = 'dashed')+
		geom_vline(xintercept = 0.1,color = 'red',linetype = 'dashed')+
		scale_color_manual(values = color_map, guide = FALSE)

	forest_plot_data = data.table(
		tissue = colnames(beta),
		beta = unlist(beta[rownames(beta)==id,]),
		se = unlist(se[rownames(se)==id,])
		)
	forest_plot_data = forest_plot_data[!is.na(beta)&!is.na(se),]

	p2 = ggplot(forest_plot_data,aes(x = tissue, y = beta, ymin = beta - 1.96*se, ymax = beta + 1.96*se, color = tissue))+
		geom_linerange(size=1.5) + 
		geom_point(size=3) + 
		ylab('Effect size') + 
		xlab('') + 
		geom_hline(yintercept = 0, color = 'red', linetype = 'dashed') + 
		theme(axis.text.y = element_text(color = ifelse(forest_plot_data$tissue == 'HCASMC','purple','black'))) + 
		scale_color_manual(values = color_map, guide = FALSE) + 
		coord_flip()

	p3 = plot_grid(p2,p1,labels=c('A','B'),align='h',rel_widths=c(2,1))
	return(p3)
}

save_tseQTL = function(tseQTL,out_fn){
	fwrite(tseQTL,out_fn,sep='\t')
}

metasoft = read_metasoft(in_file)
mvalue = subset_mvalue(metasoft)
pvalue = subset_pvalue(metasoft)

metasoft_input = read_metasoft_input(metasoft_input_fn)
beta = subset_beta(metasoft_input)
se = subset_se(metasoft_input)

tseQTL = select_tseQTL(mvalue,tissue='HCASMC',cutoff2 = 0.1)
tseQTL = subset_protein_coding_and_lncRNA(tseQTL)
color_map = get_color_map()

foreach (id = tseQTL$id)%do%{
	p = make_forestPM_plot(mvalue,pvalue,beta,se,id,color_map)
	save_plot(sprintf('%s/%s.pdf',fig_dir,id),p,base_width=8,base_height=8)
}

save_tseQTL(tseQTL,sprintf('%s/tseQTL.tsv',out_dir))


for (cutoff2 in c(0.1,0.3,0.5,0.7,0.9)){
	tseQTL = select_tseQTL(mvalue,tissue='HCASMC',cutoff2 = cutoff2)
	tseQTL = subset_protein_coding_and_lncRNA(tseQTL)
	save_tseQTL(tseQTL,sprintf('%s/tseQTL_cutoff%s.tsv',out_dir,cutoff2))
}