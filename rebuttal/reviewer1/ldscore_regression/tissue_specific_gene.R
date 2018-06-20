library(data.table)

# Variables:
gtex_rpkm_fn='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6p_All_Tissues_gene_rpkm_FOR_QC_ONLY.gct'
gtex_sample_annotation_fn='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/sample_annotations/GTEx_Analysis_2015-01-12_Annotations_SampleAttributesDS.txt'
hcasmc_rpkm_fn='../processed_data/rnaseq/preprocess/combine_rpkm/combined.rpkm'
gene_annotation_fn='../data/gtex/gencode.v19.genes.v6p.hg19.bed'
out_dir='../processed_data/rebuttal/reviewer1/ldscore_regression/tissue_specific_gene/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

# Functions: 
# Select tissue specific genes (> 4sd above the mean):
select_tissue_specific_genes_4sd=function(kept_tissue_rpkm,out_dir){
	kept_tissue_rpkm_sd=apply(kept_tissue_rpkm,1,sd)
	kept_tissue_rpkm_mean=apply(kept_tissue_rpkm,1,mean)
	for (tissue in names(kept_tissue_rpkm)){
		message(tissue)
		tissue_specific_gene=names(which(kept_tissue_rpkm[,tissue]>kept_tissue_rpkm_mean+4*kept_tissue_rpkm_sd))
		tissue_specific_gene_annotation=gene_annotation[Name%in%tissue_specific_gene,]
		fwrite(tissue_specific_gene_annotation,sprintf('%s/%s.4sd.txt',out_dir,tissue),sep='\t')
	}
}

# Select top tissue specific genes:
select_tissue_specific_genes_topN=function(kept_tissue_rpkm,out_dir,topN_list){
	kept_tissue_rpkm_sd=apply(kept_tissue_rpkm,1,sd)
	kept_tissue_rpkm_mean=apply(kept_tissue_rpkm,1,mean)
	standard_rpkm=(kept_tissue_rpkm-kept_tissue_rpkm_mean)/kept_tissue_rpkm_sd


	# Select top N tissue-specific genes:
	for (tissue in names(kept_tissue_rpkm)){
		message(tissue)
		tissue_standard_rpkm=data.table(Name=rownames(standard_rpkm),rpkm=standard_rpkm[,tissue])
		setorder(tissue_standard_rpkm,-rpkm,na.last=TRUE)

		for (n in topN_list){
			tissue_specific_gene=tissue_standard_rpkm[1:n,Name]
			tissue_specific_gene_annotation=gene_annotation[Name%in%tissue_specific_gene,]
			fwrite(tissue_specific_gene_annotation,sprintf('%s/%s.top%s.txt',out_dir,tissue,n),sep='\t')
		}
	}
}

# Read GTEx RPKMs: 
gtex=fread(gtex_rpkm_fn)
gtex_sample_annotation=fread(gtex_sample_annotation_fn)[,list(sample_id=SAMPID,tissue=SMTSD)]


# Read HCASMC RPKM:
hcasmc=fread(hcasmc_rpkm_fn,header=TRUE)


# Calculate median RPKM per tissue: 
gtex_long=melt(gtex,id.vars=c('Name','Description'),variable.name='sample_id',value.name='rpkm')
gtex_long=merge(gtex_long,gtex_sample_annotation,by='sample_id')
gtex_median_rpkm=gtex_long[,list(median_rpkm=median(rpkm)),by=c('Name','tissue')]


hcasmc_long=melt(hcasmc,id.vars=c('Name'),variable.name='sample_id',value.name='rpkm')
hcasmc_long$tissue='HCASMC'
hcasmc_median_rpkm=hcasmc_long[,list(median_rpkm=median(rpkm)),by=c('Name','tissue')]


# Combined GTEx and HCASMC RPKMs: 
combined_rpkm_long=rbind(hcasmc_median_rpkm,gtex_median_rpkm)
# fwrite(combined_rpkm_long,sprintf('%s/combined_rpkm_long.txt',out_dir),sep='\t')
# combined_rpkm_long=fread(sprintf('%s/combined_rpkm_long.txt',out_dir))

# Read gene annotation: 
gene_annotation=fread(gene_annotation_fn,select=c(5:7),col.names=c('Name','Description','type'))


# Subset to protein coding genes:
protein_coding_genes=gene_annotation[type=='protein_coding',Name]
combined_rpkm_long=combined_rpkm_long[Name%in%protein_coding_genes]


# Select independent tissue (spearman's R < 0.96):
combined_rpkm=dcast(combined_rpkm_long,Name~tissue,value.var='median_rpkm')
hcasmc_col=which(colnames(combined_rpkm)=='HCASMC')
coronary_artery_col=which(colnames(combined_rpkm)=='Artery - Coronary')
neworder=c(hcasmc_col,coronary_artery_col,seq(ncol(combined_rpkm))[-c(hcasmc_col,coronary_artery_col)])
setcolorder(combined_rpkm,neworder)


setDF(combined_rpkm)
rownames(combined_rpkm)=combined_rpkm$Name
combined_rpkm$Name=NULL

combined_rpkm_corr=cor(combined_rpkm)
diag(combined_rpkm_corr)=0 # set diagonal to 0 to keep current tissue in each iteration.
combined_rpkm_corr[3:4,3:4] # Pearson correlation = 0.9899669

threshold=0.96
n_tissue_kept=0
n_tissue_remaining=nrow(combined_rpkm_corr)
neighboring_tissue=list()
while(n_tissue_remaining>0){
	n_tissue_kept=n_tissue_kept+1

	tissue_to_remove=names(which(combined_rpkm_corr[n_tissue_kept,]>=threshold))
	neighboring_tissue[[rownames(combined_rpkm_corr)[n_tissue_kept]]]=tissue_to_remove

	tissue_to_keep=which(combined_rpkm_corr[n_tissue_kept,]<threshold)
	combined_rpkm_corr=combined_rpkm_corr[tissue_to_keep,tissue_to_keep]
	
	n_tissue_remaining=nrow(combined_rpkm_corr)-n_tissue_kept
}

tissue_kept=colnames(combined_rpkm_corr)
tissue_kept[3] = "Adipose - Visceral (Omentum)"
fwrite(data.table(tissue_kept),sprintf('%s/kept_tissue.txt',out_dir),row.names=FALSE,col.names=FALSE)
saveRDS(neighboring_tissue,sprintf('%s/neighboring_tissue.rds',out_dir))


#----------- Use all independent tissues -----------#
# Select tissue-specific genes (greater than 4 s.d.): 
kept_tissue_rpkm=combined_rpkm[tissue_kept]
out_dir_with_sm=sprintf('%s/tissue_specific_gene/',out_dir)
if(!dir.exists(out_dir_with_sm)) {dir.create(out_dir_with_sm)}
select_tissue_specific_genes_4sd(kept_tissue_rpkm,out_dir_with_sm)


# Select top tissue-specific genes:
topN_list=c(200,500,1000,2000,4000,8000)
select_tissue_specific_genes_topN(kept_tissue_rpkm,out_dir_with_sm,topN_list)


#----------- Use all independent tissues minus smooth muscle -----------#
# Remove smooth muscle tissues: 
smooth_muscle_tissue_list=c('Cervix - Endocervix','Colon - Sigmoid','Esophagus - Mucosa','Vagina','Stomach')
kept_tissue_rpkm_no_sm=kept_tissue_rpkm[,!(names(kept_tissue_rpkm)%in%smooth_muscle_tissue_list)]


# Select tissue-specific genes (greater than 4 s.d.): 
out_dir_no_sm=sprintf('%s/tissue_specific_gene_no_sm/',out_dir)
if(!dir.exists(out_dir_no_sm)) {dir.create(out_dir_no_sm)}
select_tissue_specific_genes_4sd(kept_tissue_rpkm_no_sm,out_dir_no_sm)


# Select top tissue-specific genes:
select_tissue_specific_genes_topN(kept_tissue_rpkm_no_sm,out_dir_no_sm,topN_list)


#----------- Use all independent tissues minus smooth muscle and blood -----------#
# Remove blood: 
kept_tissue_rpkm_no_sm_no_blood=kept_tissue_rpkm_no_sm[,names(kept_tissue_rpkm_no_sm)!='Whole Blood']


# Select tissue-specific genes (greater than 4 s.d.): 
out_dir_no_sm_no_blood=sprintf('%s/tissue_specific_gene_no_sm_no_blood/',out_dir)
if(!dir.exists(out_dir_no_sm_no_blood)) {dir.create(out_dir_no_sm_no_blood)}
select_tissue_specific_genes_4sd(kept_tissue_rpkm_no_sm_no_blood,out_dir_no_sm_no_blood)


# Select top tissue-specific genes:
select_tissue_specific_genes_topN(kept_tissue_rpkm_no_sm_no_blood,out_dir_no_sm_no_blood,topN_list)