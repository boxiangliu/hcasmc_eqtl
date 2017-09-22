library(data.table)

# Variables:
gtex_rpkm_fn='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6p_All_Tissues_gene_rpkm_FOR_QC_ONLY.gct'
gtex_sample_annotation_fn='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/sample_annotations/GTEx_Analysis_2015-01-12_Annotations_SampleAttributesDS.txt'
hcasmc_rpkm_fn='../processed_data/rnaseq/preprocess/combine_rpkm/combined.rpkm'
gene_annotation_fn='../data/gtex/gencode.v19.genes.v6p.hg19.bed'
out_dir='../processed_data/gwas_gene_overlap/ldscore_regression/tissue_specific_gene/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}


# Read GTEx RPKMs: 
gtex=fread(gtex_rpkm_fn)
gtex_sample_annotation=fread(gtex_sample_annotation_fn)[,
	list(sample_id=SAMPID,tissue=SMTSD)]


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


# Read gene annotation: 
gene_annotation=fread(gene_annotation_fn,select=c(5:7),col.names=c('Name','Description','type'))


# Subset to protein coding genes:
protein_coding_genes=gene_annotation[type=='protein_coding',Name]
combined_rpkm_long=combined_rpkm_long[Name%in%protein_coding_genes]


# Select independent tissue (spearman's R < 0.96):
combined_rpkm=dcast(combined_rpkm_long,Name~tissue,value.var='median_rpkm')
hcasmc_col=which(colnames(combined_rpkm)=='HCASMC')
setcolorder(combined_rpkm,c(hcasmc_col,1:(hcasmc_col-1),(hcasmc_col+1):ncol(combined_rpkm)))

setDF(combined_rpkm)
rownames(combined_rpkm)=combined_rpkm$Name
combined_rpkm$Name=NULL

combined_rpkm_corr=cor(combined_rpkm)
diag(combined_rpkm_corr)=0 # set diagonal to 0 to keep current tissue in each iteration.

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
fwrite(data.table(tissue_kept),sprintf('%s/kept_tissue.txt',out_dir),row.names=FALSE,col.names=FALSE)
saveRDS(neighboring_tissue,sprintf('%s/neighboring_tissue.rds',out_dir))


# Select tissue-specific genes: 
kept_tissue_rpkm=combined_rpkm[tissue_kept]
kept_tissue_rpkm_sd=apply(kept_tissue_rpkm,1,sd)
kept_tissue_rpkm_mean=apply(kept_tissue_rpkm,1,mean)
if(!dir.exists(sprintf('%s/tissue_specific_gene/',out_dir))) {dir.create(sprintf('%s/tissue_specific_gene/',out_dir))}

for (tissue in names(combined_rpkm)){
	tissue_specific_gene=names(which(combined_rpkm[,tissue]>kept_tissue_rpkm_mean+4*kept_tissue_rpkm_sd))
	tissue_specific_gene_rpkm=combined_rpkm[rownames(combined_rpkm)%in%tissue_specific_gene,] 
	fwrite(tissue_specific_gene_rpkm,sprintf('%s/tissue_specific_gene/%s.rpkm',out_dir,tissue),sep='\t')

	tissue_specific_gene_annotation=gene_annotation[Name%in%tissue_specific_gene,]
	fwrite(tissue_specific_gene_annotation,sprintf('%s/tissue_specific_gene/%s.txt',out_dir,tissue),sep='\t')
}
