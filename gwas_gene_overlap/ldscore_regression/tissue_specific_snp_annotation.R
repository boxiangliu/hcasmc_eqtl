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

# Read and merge PLINK bim files: 
bim_fn_list=list.files(plink_dir,'bim',full.names=TRUE)
bim=foreach(i=seq_along(bim_fn_list),.combine='rbind')%dopar%{
	bim_fn=bim_fn_list[i]
	message(bim_fn)
	fread(bim_fn,col.names=
		c('CHR','SNP','CM','BP','A1','A2'))
}


# Read gene annotation:
gene_annotation=fread(sprintf('zcat %s',gene_annotation_fn),
	select=c(1,3:5,9),
	col.names=c('CHR','TYPE','START','END','ANNOT'))
gene_annotation[,GENE_ID:=str_extract(ANNOT,'(?<=gene_id \\")(ENSG[0-9\\.]+)(?=\\";)')]
gene_annotation[,ANNOT:=NULL]


# Read tissue-specific gene: 
tissue_specific_gene_fn_list=list.files(tissue_specific_gene_dir,'txt',full.names=TRUE)

tissue_specific_gene=foreach(i=seq_along(tissue_specific_gene_fn_list),.combine='rbind')%dopar%{
	tissue_specific_gene_fn=tissue_specific_gene_fn_list[i]
	tissue=str_replace(basename(tissue_specific_gene_fn),'.txt','')
	data.table(
		fread(tissue_specific_gene_fn_list[i],select=1,col.names='GENE_ID'),
		TISSUE=tissue)
}


# Merge exon intervals and tissue specific genes:
tissue_specific_exon_interval=merge(gene_annotation[TYPE=='exon',],tissue_specific_gene,by='GENE_ID')
setorder(tissue_specific_exon_interval,TISSUE,CHR,START)
stopifnot(tissue_specific_exon_interval[,all(END>=START)])
tissue_specific_exon_interval[,c('START','END'):=list(START-1000,END+1000)]


# Merge tissue-specific exon interval and bim file: 
bim[,c('CHR','START','END'):=list(as.character(CHR),BP,BP)]
setkey(bim,CHR,START,END)
setkey(tissue_specific_exon_interval,CHR,START,END)
tissue_list=sort(unique(tissue_specific_exon_interval$TISSUE))

make_annot=function(t){
	message(t)
	annot=foverlaps(bim,tissue_specific_exon_interval[TISSUE==t])
	annot[,'ANNOT':=ifelse(is.na(GENE_ID),0.0,1.0)]
	list(tissue=t,annot=unique(annot[,list(CHR,BP,SNP,CM,ANNOT)]))
}

annot_list=mclapply(tissue_list,make_annot,mc.preschedule=FALSE)



# Output tissue-specific annotation by chromosome: 
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


# Combine tissue-specific annotations for independent tissues:
kept_tissue=fread(kept_tissue_fn,sep='\t',header=FALSE,col.names='tissue')
kept_tissue[tissue=='Artery - Aorta',tissue:="Artery - Coronary"]
kept_tissue[,tissue:=str_replace_all(tissue,' ','_')]

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


# Output merged annotation by chromosome: 
if(!dir.exists(sprintf('%s/merged',out_dir))){dir.create(sprintf('%s/merged',out_dir))}

foreach(chr=seq(22))%dopar%{
	message(chr)
	out_fn=gzfile(sprintf('%s/merged/merged.%s.annot.gz',out_dir,chr))
	write.table(merged_annot[CHR==chr],out_fn,sep='\t',quote=FALSE,col.names=T,row.names=F)
}


