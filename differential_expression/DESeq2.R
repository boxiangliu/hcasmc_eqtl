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


# Function:
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


# Read HCASMC read count:
hcasmc_read_count=fread(hcasmc_read_count_fn,header=TRUE)


# Subset to protein coding genes and lncRNAs:
gencode=fread('/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/gencode.v19.genes.v6p.hg19.bed',select=c(5,6,7),col.names=c('gene_id','gene_name','type'))
coding_and_lncRNA=gencode[type%in%c('lincRNA','protein_coding'),gene_id]
hcasmc_read_count=hcasmc_read_count[Name%in%coding_and_lncRNA]
hcasmc_read_count_mat=as.matrix(hcasmc_read_count[,3:ncol(hcasmc_read_count)])
rownames(hcasmc_read_count_mat)=hcasmc_read_count[,Name]


# Get HCASMC sample list: 
hcasmc_subject_list=colnames(hcasmc_read_count)[3:ncol(hcasmc_read_count)]


# Read HCASMC covariates: 
hcasmc_covariate=readWorksheet(loadWorkbook(hcasmc_covariate_fn),sheet=1)
setDT(hcasmc_covariate)
hcasmc_covariate=unique(hcasmc_covariate[DNA.New.Name%in%hcasmc_subject_list,list(SUBJID=DNA.New.Name,sex=Sex,age=Age,race=Genomic_Ethnicity)])


# Read GTEx covariates:
gtex_covariate=fread(gtex_covariate_fn)[,list(SUBJID,age=AGE,sex=GENDER,race=RACE,ethnicity=ETHNCTY)]
race_numeric_to_letter=c('Asian','AA','Caucasian','Indian','Unknown')
gtex_covariate[,race:=race_numeric_to_letter[race]]
gtex_covariate[,race:=ifelse(ethnicity==1,'Hispanic',race)]
gtex_covariate[,ethnicity:=NULL]
sex_numeric_to_letter=c('M','F')
gtex_covariate[,sex:=sex_numeric_to_letter[sex]]


# Tissue list: 
tissue_list=c("Artery - Aorta","Artery - Coronary","Artery - Tibial", 'Heart - Atrial Appendage', 'Heart - Left Ventricle','Cells - Transformed fibroblasts')
tissue_to_id_map=fread('/srv/persistent/bliu2/HCASMC_eQTL/data//gtex/gtex_tissue_colors.txt',select=c(1,3))


# Iterate across tissues:
container=list()
for (i in 1:length(tissue_list)){
	tissue=tissue_list[i]


	# Convert tissue name to tissue ID:
	tissue_id=tissue_to_id_map[tissue_site_detail==tissue,tissue_site_detail_id]
	message('INFO - ', tissue_id)


	# Read GTEx read count: 
	message('INFO - reading GTEx read counts...')
	gtex_read_count_fn=sprintf('%s/%s.reads.txt',gtex_read_count_dir,tissue_id)
	gtex_read_count=fread(gtex_read_count_fn)
	gtex_read_count=gtex_read_count[Gene%in%coding_and_lncRNA]


	# Get list of subject in the subsample:
	subject_fn=sprintf('%s/%s.txt',subject_dir,tissue)
	subject_list=fread(subject_fn,header=FALSE,col.names='subject')


	# Get read count for subject in subsample:
	gtex_read_count_downsample_mat=as.matrix(gtex_read_count[,unlist(subject_list),with=F])
	rownames(gtex_read_count_downsample_mat)=gtex_read_count[,Gene]


	# Read GTEx covariate:
	message('INFO - getting GTEx covariates...')
	tissue_covariate=gtex_covariate[SUBJID%in%unlist(subject_list)]


	# Merge count table:
	message('INFO - merging HCASMC and GTEx count table...')
	row_idx=match(rownames(hcasmc_read_count_mat),rownames(gtex_read_count_downsample_mat))
	gtex_read_count_downsample_mat_reordered=gtex_read_count_downsample_mat[row_idx,]
	stopifnot(rownames(hcasmc_read_count_mat)==rownames(gtex_read_count_downsample_mat_reordered))
	read_count_mat=cbind(hcasmc_read_count_mat,gtex_read_count_downsample_mat_reordered)


	# Merge covariates: 
	message('INFO - merging HCASCM and GTEx covariates...')
	hcasmc_covariate$tissue='HCASMC'
	tissue_covariate$tissue=tissue
	covariate=rbind(hcasmc_covariate,tissue_covariate)


	# Create DESeq dataset: 
	dds=DESeqDataSetFromMatrix(countData=read_count_mat,colData=covariate,design = ~ tissue+sex+race)


	# Perform DESeq2:
	message('INFO - running DESeq2')
	dds=DESeq(dds,parallel=TRUE)
	container[[tissue]]=dds
}


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
saveRDS(meta_result,sprintf('%s/meta_result.rds',out_dir))


# Count the number of DE genes:
fdr_threshold=c(0.001,0.01,0.05)
n_sig=foreach(t=names(meta_result),.combine='rbind')%do%{
	foreach(ft=fdr_threshold,.combine='rbind')%do%{
		data.table(tissue=t,FDR=ft,n=sum(meta_result[[t]]$padj<ft,na.rm=TRUE))
	}

}
fwrite(n_sig,sprintf('%s/n_sig.txt',out_dir),sep='\t')


# Plot the number of significant DE genes:
n_sig[,FDR:=as.character(FDR)]
n_sig[,tissue:=factor(tissue,c('Fibroblast','Heart','Artery'))]
p1=ggplot(n_sig,aes(tissue,n,fill=FDR))+geom_bar(stat='identity',position=position_dodge())+ylab('DE genes')+xlab('')
p2=ggplot(n_sig,aes(tissue,n,fill=FDR,label=n))+geom_bar(stat='identity',position=position_dodge())+ylab('DE genes')+xlab('')+geom_text(position=position_dodge(0.9),vjust= -0.5)
pdf(sprintf('%s/n_sig.pdf',fig_dir))
p1;p2
dev.off()

