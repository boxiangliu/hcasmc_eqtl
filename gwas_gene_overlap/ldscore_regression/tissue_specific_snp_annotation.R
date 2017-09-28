library(data.table)
library(foreach)
library(doMC)
registerDoMC(20)
library(stringr)

# Variables:
tissue_specific_gene_dir='../processed_data/gwas_gene_overlap/ldscore_regression/tissue_specific_gene/tissue_specific_gene/'
kept_tissue_fn='../processed_data/gwas_gene_overlap/ldscore_regression/tissue_specific_gene/kept_tissue.txt'
gene_annotation_fn='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/reference_files/gencode.v19.genes.v6p.patched_contigs.gtf.gz'
plink_dir='/srv/persistent/bliu2/shared/ldscore/1000G_plinkfiles/'
out_dir='../processed_data/gwas_gene_overlap/ldscore_regression/tissue_specific_snp_annotation/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

# Functions:
read_tissue_specific_gene=function(tissue_specific_gene_fn_list){
	tissue_specific_gene=foreach(i=seq_along(tissue_specific_gene_fn_list),.combine='rbind')%dopar%{
		tissue_specific_gene_fn=tissue_specific_gene_fn_list[i]
		tissue=str_split_fixed(basename(tissue_specific_gene_fn),'\\.',2)[,1]
		data.table(
			fread(tissue_specific_gene_fn_list[i],select=1,col.names='GENE_ID'),
			TISSUE=tissue)
	}
	return(tissue_specific_gene)
}

merge_exon_and_tissue_specific_gene=function(gene_annotation,tissue_specific_gene){
	tissue_specific_exon_interval=merge(gene_annotation[TYPE=='exon',],tissue_specific_gene,by='GENE_ID',allow.cartesian=TRUE)
	setorder(tissue_specific_exon_interval,TISSUE,CHR,START)
	stopifnot(tissue_specific_exon_interval[,all(END>=START)])
	tissue_specific_exon_interval[,c('START','END'):=list(START-1000,END+1000)]
	return(tissue_specific_exon_interval)
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

	return(annot)
}

output_tissue_specific_annotation=function(annot_list,out_dir){
	if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
	for (i in seq_along(annot_list)){
		tissue=annot_list[[i]][['tissue']]
		message(tissue)
		annot=annot_list[[i]][['annot']]

		foreach(chr=seq(22))%dopar%{
			message(chr)
			out_fn=gzfile(sprintf('%s/%s.%s.annot.gz',out_dir,tissue,chr))
			write.table(annot[CHR==chr],out_fn,sep='\t',quote=FALSE,col.names=T,row.names=F)
		}
	}
}


merge_tissue_specific_annotations=function(annot_list){
	merged_annot=annot_list[[1]][['annot']]
	setnames(merged_annot,'ANNOT',str_replace_all(annot_list[[1]][['tissue']],' ','_'))
	for(i in 2:length(annot_list)){
		
		tissue=str_replace_all(annot_list[[i]][['tissue']],' ','_')
		
		if (tissue %in% kept_tissue$tissue){
		
			message(tissue)
		
			annot=annot_list[[i]][['annot']]

			if ('ANNOT' %in% colnames(annot)){
				setnames(annot,'ANNOT',tissue)
			}
		
			merged_annot=merge(merged_annot,annot,by=c('CHR','BP','SNP','CM'))
		}
	}
	return(merged_annot)
}


output_merged_annotations=function(merged_annot,out_dir){
	if(!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
	foreach(chr=seq(22))%dopar%{
		message(chr)
		out_fn=gzfile(sprintf('%s/merged.%s.annot.gz',out_dir,chr))
		write.table(merged_annot[CHR==chr],out_fn,sep='\t',quote=FALSE,col.names=T,row.names=F,na='0')
	}
}


driver=function(pattern,out_dir_extension){
	# Read tissue-specific gene: 
	message('INFO - Reading tissue-specific gene')
	tissue_specific_gene_fn_list=list.files(tissue_specific_gene_dir,pattern,full.names=TRUE)
	tissue_specific_gene=read_tissue_specific_gene(tissue_specific_gene_fn_list)


	# Merge exon intervals and tissue specific genes:
	message('INFO - Merging exon interval and tissue specific genes')
	tissue_specific_exon_interval=merge_exon_and_tissue_specific_gene(gene_annotation,tissue_specific_gene)
	tissue_specific_exon_interval_wide=dcast(tissue_specific_exon_interval,CHR+START+END~TISSUE,value.var='GENE_ID',fun.aggregate=length)

	# Merge tissue-specific exon interval and bim file:
	message('INFO - Merging tissue-specific exon interval and bim file')
	annot=merge_tissue_specific_exon_and_bim(tissue_specific_exon_interval_wide,bim,tissue_list)


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
bim[,c('CHR','START','END'):=list(as.character(CHR),BP,BP)]
setkey(bim,CHR,START,END)


# Read gene annotation:
gene_annotation=fread(sprintf('zcat %s',gene_annotation_fn),
	select=c(1,3:5,9),
	col.names=c('CHR','TYPE','START','END','ANNOT'))
gene_annotation[,GENE_ID:=str_extract(ANNOT,'(?<=gene_id \\")(ENSG[0-9\\.]+)(?=\\";)')]
gene_annotation[,ANNOT:=NULL]


# Read kept tissues:
kept_tissue=fread(kept_tissue_fn,sep='\t',header=FALSE,col.names='tissue')
kept_tissue[,tissue:=str_replace_all(tissue,' ','_')]


# Do the magic:
driver('4sd.txt','4sd')
driver('top200.txt','top200')
driver('top500.txt','top500')
driver('top1000.txt','top1000')
