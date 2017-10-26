library(data.table)
library(foreach)
library(doMC)
registerDoMC(15)
library(stringr)

# Variables:
args=commandArgs(T)
if (length(args)==0) {
	bed_dir='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_released_all_cell_type/'
	out_dir='../processed_data/gwas_atacseq_overlap/ldscore_regression/tissue_specific_snp_annotation/'
} else if (length(args)==2){
	bed_dir=args[1]
	out_dir=args[2]
} else {
	stop('Arg 1: bed_dir; Arg 2: out_dir')
}

plink_dir='/srv/persistent/bliu2/shared/ldscore/1000G_plinkfiles/'

if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

# Functions:
read_bed=function(bed_fn_list,kept_tissue){
	bed=foreach(i=seq_along(bed_fn_list),.combine='rbind')%dopar%{
		bed_fn=bed_fn_list[i]
		tissue=str_split_fixed(basename(bed_fn),'\\.',2)[,1]
		message(tissue)

		if (tissue %in% kept_tissue | is.null(kept_tissue)){
			data.table(
				fread(bed_fn,select=1:3,col.names=c('CHR','START','END')),
				TISSUE=tissue)
		} else {
			data.table()
		}

	}
	return(bed)
}


merge_tissue_specific_exon_and_bim=function(tissue_specific_exon_interval_wide,bim,tissue_list){
	setkey(tissue_specific_exon_interval_wide,CHR,START,END)
	temp=foverlaps(bim,tissue_specific_exon_interval_wide)
	annot=unique(temp[,c('CHR','BP','SNP','CM',tissue_list),with=F])
	
	setDF(annot)
	annot[is.na(annot)]=0
	
	setDT(annot)
	dup=which(duplicated(annot[,list(CHR,BP,SNP,CM)]))

	for (tissue in tissue_list){
		setnames(annot,tissue,'temp')
		annot[c(dup,dup-1),temp:=max(temp,na.rm=TRUE),by=c('CHR','BP','SNP','CM')]
		setnames(annot,'temp',tissue)
	}

	annot=unique(annot)
	stopifnot(nrow(annot)==nrow(annot))
	stopifnot(annot$SNP==bim$SNP)

	annot
	return(annot)
}

output_merged_annotations=function(merged_annot,out_dir){
	if(!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
	foreach(chr=seq(22))%dopar%{
		message(chr)
		out_fn=gzfile(sprintf('%s/merged.%s.annot.gz',out_dir,chr))
		write.table(merged_annot[CHR==chr],out_fn,sep='\t',quote=FALSE,col.names=T,row.names=F,na='0')
	}
}


jaccard=function(a,b){
	intersection=a&b
	union=a|b
	sum(intersection)/sum(union)
}


select_independent_tissue=function(cor_mat,threshold){
	diag(cor_mat)=0 # set diagonal to 0 to keep current tissue in each iteration.

	n_tissue_kept=0
	n_tissue_remaining=nrow(cor_mat)
	neighboring_tissue=list()
	while(n_tissue_remaining>0){
		n_tissue_kept=n_tissue_kept+1

		tissue_to_remove=names(which(cor_mat[n_tissue_kept,]>=threshold))
		neighboring_tissue[[rownames(cor_mat)[n_tissue_kept]]]=tissue_to_remove

		tissue_to_keep=which(cor_mat[n_tissue_kept,]<threshold)
		cor_mat=cor_mat[tissue_to_keep,tissue_to_keep]
		
		n_tissue_remaining=nrow(cor_mat)-n_tissue_kept
	}
	tissue_kept=colnames(cor_mat)
	return(list(tissue_kept=tissue_kept,neighboring_tissue=neighboring_tissue))
}


driver=function(in_dir,pattern,out_dir_extension){
	# Read tissue-specific gene: 
	message('INFO - Reading tissue-specific gene')
	tissue_specific_gene_fn_list=list.files(in_dir,pattern,full.names=TRUE)
	tissue_specific_gene=read_tissue_specific_gene(tissue_specific_gene_fn_list,kept_tissue$tissue)


	# Merge exon intervals and tissue specific genes:
	message('INFO - Merging exon interval and tissue specific genes')
	tissue_specific_exon_interval=merge_exon_and_tissue_specific_gene(gene_annotation,tissue_specific_gene)
	tissue_specific_exon_interval_wide=dcast(tissue_specific_exon_interval,CHR+START+END~TISSUE,value.var='GENE_ID',fun.aggregate=length)

	# Merge tissue-specific exon interval and bim file:
	message('INFO - Merging tissue-specific exon interval and bim file')
	tissue_list=sort(unique(tissue_specific_exon_interval$TISSUE))
	annot=merge_tissue_specific_exon_and_bim(tissue_specific_exon_interval_wide,bim,tissue_list)
	setnames(annot,names(annot),str_replace_all(names(annot),' ','_'))

	# Output merged annotation by chromosome:
	message('INFO - Outputing merged annotation by chromosome')
	output_merged_annotations(annot,sprintf('%s/%s/',out_dir,out_dir_extension))
}

# Read and merge PLINK bim files: 
bim_fn_list=list.files(plink_dir,'bim',full.names=TRUE)
bim=foreach(i=seq_along(bim_fn_list),.combine='rbind')%dopar%{
	bim_fn=bim_fn_list[i]
	message(bim_fn)
	fread(bim_fn,col.names=
		c('CHR','SNP','CM','BP','A1','A2'))
}
bim[,c('CHR','START','END'):=list(paste0('chr',CHR),BP,BP)]
setkey(bim,CHR,START,END)


# Read bed files: 
message('INFO - Reading tissue-specific gene')
bed_fn_list=list.files(bed_dir,pattern='bed',full.names=TRUE)
bed=read_bed(bed_fn_list,NULL)
tissue_list=unique(bed$TISSUE)
annot=foreach(tissue=tissue_list,.combine='cbind')%dopar%{
	message(tissue)
	tissue_bed=bed[TISSUE==tissue]
	setkey(tissue_bed,CHR,START,END)
	overlap=foverlaps(bim,tissue_bed)
	overlap[,temp:=!is.na(TISSUE)]
	tissue_annot=overlap[,max(temp),by=c('CHR','BP','SNP','CM')]
	stopifnot(nrow(tissue_annot)==nrow(bim)) # quick check
	setnames(tissue_annot,'V1',tissue)
	tissue_annot[,tissue,with=F]
}

temp=cbind(bim[,list(CHR,BP,SNP,CM)],annot)
temp[,CHR:=str_replace(CHR,'chr','')]
annot=temp
setnames(annot,str_replace_all(colnames(annot),' ','_'))


# Output merged annotation by chromosome:
message('INFO - Outputing merged annotation by chromosome')
output_merged_annotations(annot,sprintf('%s/all_tissue/',out_dir))

# annot=foreach(i=1:22,.combine='rbind')%dopar%{
# 	message(i)
# 	fread(sprintf('zcat ../processed_data/gwas_atacseq_overlap/ldscore_regression_2305/tissue_specific_snp_annotation/all_tissue/merged.%s.annot.gz',i))}

# Calculate jaccard similarity:
size=ncol(annot)-4
jaccard_similarity=matrix(-1,nrow=size,ncol=size)
colnames(jaccard_similarity)=rownames(jaccard_similarity)=colnames(annot)[5:ncol(annot)]
registerDoMC(48)
for (i in 5:ncol(annot)){
	temp=foreach(j=i:ncol(annot),.combine='c')%dopar%{
		message(i,',',j)
		jaccard(annot[,i,with=FALSE],annot[,j,with=FALSE])
	}
	jaccard_similarity[i-4,i:ncol(annot)-4]=temp
}
jaccard_similarity[lower.tri(jaccard_similarity)]=t(jaccard_similarity)[lower.tri(jaccard_similarity)]
out=data.table(jaccard_similarity,keep.rownames=TRUE)
setnames(out,'rn','tissue')
fwrite(out,sprintf('%s/jaccard_similarity.tsv',out_dir),sep='\t')

# out_dir='../processed_data/gwas_atacseq_overlap/ldscore_regression_2305/tissue_specific_snp_annotation/'
# out=fread(sprintf('%s/jaccard_similarity.tsv',out_dir))
# jaccard_similarity=out[,2:ncol(out)]
# setDF(jaccard_similarity)
# rownames(jaccard_similarity)=out[,tissue]

# Select indepedent tissues: 
new_order=c('HCASMC',colnames(jaccard_similarity)[colnames(jaccard_similarity)!='HCASMC'])
jaccard_similarity_hcasmc_1st=jaccard_similarity[new_order,new_order]

for (threshold in c(0.3,0.4,0.5)){
	temp=select_independent_tissue(jaccard_similarity_hcasmc_1st,threshold)

	tissue_kept=temp[[1]]
	write.table(tissue_kept,sprintf('%s/tissue_kept.jaccard_similarity_%s.txt',out_dir,threshold))

	neighboring_tissue=temp[[2]]
	saveRDS(neighboring_tissue,sprintf('%s/neighboring_tissue.jaccard_similarity_%s.rds',out_dir,threshold))

	annot_threshold=annot[,c('CHR','BP','SNP','CM',tissue_kept),with=F]
	message(sprintf('INFO - %s tissue kept',ncol(annot_threshold)-4))

	# Output independent tissue annotation by chromosom:
	output_merged_annotations(annot_threshold,sprintf('%s/jaccard_similarity_%s/',out_dir,threshold))
}