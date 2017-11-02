library(data.table)
library(DESeq2)
library(XLConnect)
library(stringr)
library(BiocParallel)
register(MulticoreParam(4))
library(foreach)
library(doMC)
registerDoMC(10)
library(metafor)
library(cowplot)
library(sva)

# Variables:
subject_dir='../processed_data/differential_expression/sample_individuals/'
gtex_read_count_dir='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6p_read_counts/'
gtex_covariate_fn='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/sample_annotations/GTEx_Analysis_2015-01-12_Annotations_SubjectPhenotypesDS.txt'
hcasmc_read_count_fn='../data/rnaseq2/read_count/rnaseqc/rnaseqc.hcasmc_eqtl.reads.gct'
hcasmc_covariate_fn='/srv/persistent/bliu2/HCASMC_eQTL/data/sample_info/sample_info.xlsx'
out_dir='../processed_data/differential_expression/DESeq2/'
fig_dir='../figures/differential_expression/DESeq2/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}


tissue_list=c("Artery - Aorta","Artery - Coronary","Artery - Tibial", 'Heart - Atrial Appendage', 'Heart - Left Ventricle','Cells - Transformed fibroblasts')
tissue_to_id_map=fread('/srv/persistent/bliu2/HCASMC_eQTL/data//gtex/gtex_tissue_colors.txt',select=c(1,3))
gencode=fread('/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/gencode.v19.genes.v6p.hg19.bed',select=c(5,6,7),col.names=c('gene_id','gene_name','type'))
coding_and_lncRNA=gencode[type%in%c('lincRNA','protein_coding'),gene_id]
gtex_covariate=read_gtex_covariate(gtex_covariate_fn)


# Function:
read_hcasmc_covariate=function(hcasmc_covariate_fn){
	hcasmc_covariate=readWorksheet(loadWorkbook(hcasmc_covariate_fn),sheet=5)
	setDT(hcasmc_covariate)
	hcasmc_covariate=unique(hcasmc_covariate[DNA%in%hcasmc_subject_list,list(SUBJID=DNA,sex=Sex,age=Age_Imputed,race=Genomic_Ethnicity)])
	return(hcasmc_covariate)
}


read_gtex_covariate=function(gtex_covariate_fn){
	gtex_covariate=fread(gtex_covariate_fn)[,list(SUBJID,age=AGE,sex=GENDER,race=RACE,ethnicity=ETHNCTY)]
	race_numeric_to_letter=c('Asian','AA','Caucasian','Indian','Unknown')
	gtex_covariate[,race:=race_numeric_to_letter[race]]
	gtex_covariate[,race:=ifelse(ethnicity==1,'Hispanic',race)]
	gtex_covariate[,ethnicity:=NULL]
	sex_numeric_to_letter=c('M','F')
	gtex_covariate[,sex:=sex_numeric_to_letter[sex]]
	return(gtex_covariate)
}


read_gtex_read_count=function(gtex_read_count_dir,tissue_id,subset=coding_and_lncRNA){
	gtex_read_count_fn=sprintf('%s/%s.reads.txt',gtex_read_count_dir,tissue_id)
	gtex_read_count=fread(gtex_read_count_fn)
	gtex_read_count=gtex_read_count[Gene%in%subset]
	return(gtex_read_count)
}


get_subject_list=function(suject_dir,tissue){
	subject_fn=sprintf('%s/%s.txt',subject_dir,tissue)
	subject_list=fread(subject_fn,header=FALSE,col.names='subject')
	return(subject_list)
}


merge_count_table=function(hcasmc_read_count_mat,gtex_read_count_downsample_mat){
	row_idx=match(rownames(hcasmc_read_count_mat),rownames(gtex_read_count_downsample_mat))
	gtex_read_count_downsample_mat_reordered=gtex_read_count_downsample_mat[row_idx,]
	stopifnot(rownames(hcasmc_read_count_mat)==rownames(gtex_read_count_downsample_mat_reordered))
	read_count_mat=cbind(hcasmc_read_count_mat,gtex_read_count_downsample_mat_reordered)
	return(read_count_mat)
}


merge_covariates=function(hcasmc_covariate,tissue_covariate,tissue1='HCASMC',tissue2=tissue){
	hcasmc_covariate$tissue=tissue1
	tissue_covariate$tissue=tissue2
	covariate=rbind(hcasmc_covariate,tissue_covariate)
	return(covariate)
}


extract_sv=function(read_count_mat,covariate){
	filter = apply(read_count_mat, 1, function(x) length(x[x>5])>=2)
	filtered = read_count_mat[filter,]
	mod1 = model.matrix(~sex+race+tissue,covariate)
	mod0 = cbind(mod1[,1])
	svseq=svaseq(filtered,mod1,mod0)$sv
	colnames(svseq)=paste0('SV',1:ncol(svseq))
	return(svseq)
}


meta_analyze=function(res,tissue_group){
	meta_res=foreach(gid=unique(res[,gene_id]),.combine='rbind')%dopar%{
		message('INFO - ', gid)
		dat=res[gene_id==gid & tissue%in%tissue_group,list(log2FoldChange,lfcSE)]
		if (all(is.na(dat$log2FoldChange))){
			message('INFO - all values are missing')
			data.table()
		} else {
			temp=rma(log2FoldChange, lfcSE, data=dat, method="REML")
			data.table(gene_id=gid,beta=temp$beta[,1],se=temp$se,pval=temp$pval,het_p=temp$QEp)
		}
	}
	meta_res[,padj:=p.adjust(pval,method='BH')]
	meta_res
}


count_de_genes=function(result,fdr_threshold=c(0.001,0.01,0.05)){
	# Note that result is a list.
	n_sig=foreach(t=names(result),.combine='rbind')%do%{
		foreach(ft=fdr_threshold,.combine='rbind')%do%{
			data.table(tissue=t,FDR=ft,n=sum(result[[t]]$padj<ft,na.rm=TRUE))
		}
	}
	return(n_sig)
}


reorder_covariate_by_read_count_matrix=function(covariate,read_count){
	covariate[match(colnames(read_count),covariate$SUBJID)]
}


compare_two_tissues=function(tissue1,tissue2){

	# Read read count and covariate:
	temp=foreach(tissue=c(tissue1,tissue2))%do%{
		tissue_id=tissue_to_id_map[tissue_site_detail==tissue,tissue_site_detail_id]
		gtex_read_count=read_gtex_read_count(gtex_read_count_dir,tissue_id)


		# Get list of subject in the subsample:
		subject_list=get_subject_list(subject_dir,tissue)


		# Get read count for subject in subsample:
		gtex_read_count_downsample_mat=as.matrix(gtex_read_count[,unlist(subject_list),with=F])
		rownames(gtex_read_count_downsample_mat)=gtex_read_count[,Gene]


		# Get covariate:
		tissue_covariate=gtex_covariate[SUBJID%in%unlist(subject_list)]
		tissue_covariate=reorder_covariate_by_read_count_matrix(tissue_covariate,gtex_read_count_downsample_mat)
		list(gtex_read_count_downsample_mat,tissue_covariate)
	}


	# Merge count table:
	message('INFO - merging HCASMC and GTEx count table...')
	read_count_mat=merge_count_table(temp[[1]][[1]],temp[[2]][[1]])


	# Merge covariates: 
	message('INFO - merging HCASCM and GTEx covariates...')
	covariate=merge_covariates(temp[[1]][[2]],temp[[2]][[2]],tissue1,tissue2)


	# Create DESeq dataset: 
	dds=DESeqDataSetFromMatrix(countData=read_count_mat,colData=covariate,design = ~ tissue+sex+race)


	# Perform DESeq2:
	message('INFO - running DESeq2')
	dds=DESeq(dds,parallel=TRUE)
	return(dds)
}


# Read HCASMC read count:
hcasmc_read_count=fread(hcasmc_read_count_fn,header=TRUE)


# Subset to protein coding genes and lncRNAs:
hcasmc_read_count=hcasmc_read_count[Name%in%coding_and_lncRNA]
hcasmc_read_count_mat=as.matrix(hcasmc_read_count[,3:ncol(hcasmc_read_count)])
rownames(hcasmc_read_count_mat)=hcasmc_read_count[,Name]


# Get HCASMC sample list: 
hcasmc_subject_list=colnames(hcasmc_read_count)[3:ncol(hcasmc_read_count)]


# Read HCASMC covariates:
hcasmc_covariate=read_hcasmc_covariate(hcasmc_covariate_fn)


# Read GTEx covariates:
gtex_covariate=read_gtex_covariate(gtex_covariate_fn)


# Iterate across tissues:
container=list()
for (i in 1:length(tissue_list)){
	tissue=tissue_list[i]


	# Convert tissue name to tissue ID:
	tissue_id=tissue_to_id_map[tissue_site_detail==tissue,tissue_site_detail_id]
	message('INFO - ', tissue_id)


	# Read GTEx read count: 
	message('INFO - reading GTEx read counts...')
	gtex_read_count=read_gtex_read_count(gtex_read_count_dir,tissue_id)


	# Get list of subject in the subsample:
	subject_list=get_subject_list(subject_dir,tissue)


	# Get read count for subject in subsample:
	gtex_read_count_downsample_mat=as.matrix(gtex_read_count[,unlist(subject_list),with=F])
	rownames(gtex_read_count_downsample_mat)=gtex_read_count[,Gene]


	# Read GTEx covariate:
	message('INFO - getting GTEx covariates...')
	tissue_covariate=gtex_covariate[SUBJID%in%unlist(subject_list)]


	# Merge count table:
	message('INFO - merging HCASMC and GTEx count table...')
	read_count_mat=merge_count_table(hcasmc_read_count_mat,gtex_read_count_downsample_mat)

	# Merge covariates: 
	message('INFO - merging HCASCM and GTEx covariates...')
	covariate=merge_covariates(hcasmc_covariate,tissue_covariate)
	covariate=reorder_covariate_by_read_count_matrix(covariate,read_count_mat)
	stopifnot(covariate$SUBJID==colnames(read_count_mat))


	# Extract hidden confounders:
	svseq=extract_sv(read_count_mat,covariate)
	covariate=cbind(covariate,svseq)


	# Create DESeq dataset:
	form=as.formula(paste0('~sex+race+tissue',paste0('+SV',1:ncol(svseq),collapse='')))
	dds=DESeqDataSetFromMatrix(countData=read_count_mat,colData=covariate,design = form)


	# Perform DESeq2:
	message('INFO - running DESeq2')
	dds=DESeq(dds,parallel=TRUE)
	container[[tissue]]=dds
}
saveRDS(container,sprintf('%s/dds.rds',out_dir))

# Get results for each tissue:
res=foreach(i=1:length(tissue_list),.combine='rbind')%dopar%{
	message('INFO - ', tissue_list[i])
	res=results(container[[i]],contrast=c('tissue','HCASMC',tissue_list[i]))
	temp=as.data.table(as.data.frame(res),keep.rownames=TRUE)
	temp[,tissue:=tissue_list[i]]
	temp
}
setnames(res,'rn','gene_id')
res=res[baseMean>0] # remove genes with zero reads.


# Meta-analysis with random effect:
meta_artery=meta_analyze(res,c("Artery - Aorta","Artery - Coronary","Artery - Tibial"))
meta_heart=meta_analyze(res,c('Heart - Atrial Appendage', 'Heart - Left Ventricle'))
meta_fibroblast=res[tissue=='Cells - Transformed fibroblasts',list(gene_id,beta=log2FoldChange,se=lfcSE,pval=pvalue)]
meta_fibroblast[,padj:=p.adjust(pval,method='BH')]
meta_result=list(Artery=meta_artery,Heart=meta_heart,Fibroblast=meta_fibroblast)


# Save result to output:
saveRDS(meta_result,sprintf('%s/meta_result.rds',out_dir))


# Convert result to GSEA RNK format:
gene_annotation_fn='/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/gencode.v19.genes.v6p.hg19.bed'
gene_annotation=fread(gene_annotation_fn,select=5:6,col.names=c('gene_id','gene_name'))
gene_annotation_vect=gene_annotation$gene_name
names(gene_annotation_vect)=gene_annotation$gene_id


for (i in list('meta_artery','meta_heart','meta_fibroblast')){
	dt=get(i)
	dt[,gene_name:=gene_annotation_vect[gene_id]]
	dt[,logp:=-log10(pval)]
	dt[,z:=beta/se]
	fwrite(dt[,list(gene_name,z)],sprintf('%s/%s.rnk',out_dir,i),sep='\t',col.names=FALSE)
	fwrite(dt,sprintf('%s/%s.txt',out_dir,i),sep='\t')
}


# Count the number of DE genes:
n_sig=count_de_genes(meta_result)
fwrite(n_sig,sprintf('%s/n_sig.txt',out_dir),sep='\t')


# Plot the number of significant DE genes:
n_sig[,FDR:=as.character(FDR)]
n_sig[,tissue:=factor(tissue,c('Fibroblast','Heart','Artery'))]
p1=ggplot(n_sig,aes(tissue,n,fill=FDR))+geom_bar(stat='identity',position=position_dodge())+ylab('DE genes')+xlab('')
p2=ggplot(n_sig,aes(tissue,n,fill=FDR,label=n))+geom_bar(stat='identity',position=position_dodge())+ylab('DE genes')+xlab('')+geom_text(position=position_dodge(0.9),vjust= -0.5)
pdf(sprintf('%s/n_sig.pdf',fig_dir))
p1;p2
dev.off()


# # Differential expression between fibroblast and heart,
# # for quality control:
# res_fibroblast_vs_atrial_appendage=results(compare_two_tissues('Cells - Transformed fibroblasts','Heart - Atrial Appendage'))
# res_fibroblast_vs_coronary_artery=results(compare_two_tissues('Cells - Transformed fibroblasts','Artery - Coronary'))

# # Count the number of DE genes:
# n_sig=count_de_genes(list(fibroblast_vs_atrial_appendage=res_fibroblast_vs_coronary_artery,fibroblast_vs_coronary_artery=res_fibroblast_vs_coronary_artery))




