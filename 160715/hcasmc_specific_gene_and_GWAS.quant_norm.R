library(data.table)
library(dtplyr)
library(dplyr)
library(stringr)
library(cowplot)

# Function:
reformat=function(x){
	if (class(x$CHR_POS)!='numeric') {x[,CHR_POS:=as.integer(CHR_POS)]}
	y=x[!is.na(CHR_POS),]
	y=y[CHR_ID!='X',]
	return(y)
}

liftOver=function(x,chain='/srv/persistent/bliu2/shared/chain_files/hg38ToHg19.over.chain.gz'){
	x[,rowid:=1:nrow(x)]
	bedold=x[,.(CHR_ID,CHR_POS-1,CHR_POS,rowid)]
	bedold[,CHR_ID:=paste0('chr',CHR_ID)]
	dir.create('/srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_gene/')
	bedold_file='/srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_gene/gwas_catalog.hg38.bed'
	bednew_file='/srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_gene/gwas_catalog.hg19.bed'
	fwrite(bedold,bedold_file,sep='\t',quote=F,col.names=F)
	system(paste('/srv/persistent/bliu2/tools/ucsc_tools/liftOver',bedold_file,chain,bednew_file,'/srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_gene/gwas_catalog.unmapped.bed'))
	bednew=fread(bednew_file)[,c(3,4)]
	setnames(bednew,c('CHR_POS_19','rowid'))
	y=merge(bednew,x,by='rowid')
	y[,rowid:=NULL]
	unlink('/srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_gene/',recursive=T)
	return(y)
}


prune=function(x,min_dist=1e6){
	stopifnot(class(x$CHR_POS)=='integer')
	stopifnot(class(x$CHR_ID)=='character')
	studies=unique(x[,.(PUBMEDID,MAPPED_TRAIT)])
	z=data.table()
	for (i in 1:nrow(studies)){
		pubmedid=studies[i,PUBMEDID]
		mapped_trait=studies[i,MAPPED_TRAIT]
		y=x[ (PUBMEDID==pubmedid) & (MAPPED_TRAIT==mapped_trait),]
		y[,rank:=rank(`P-VALUE`)]
		if (nrow(y)==1){
			z=rbind(z,y)
		} else {
			setkey(y,rank)
			keep=rep(TRUE,nrow(y))
			for (j in 2:nrow(y)){
				for (k in 1:(j-1)){
					c1=y[j,CHR_ID]==y[k,CHR_ID]
					c2=(abs(y[j,CHR_POS]-y[k,CHR_POS])<min_dist)
					if (c1&c2){
						keep[j]=FALSE
					}
				}
			}
			z=rbind(z,y[keep,])
		}
	}
	return(z)
}


rmdupvar=function(x){
	y=x[,.(CHR_ID,CHR_POS_19,`P-VALUE`)]
	y=as.data.table(y%>%group_by(CHR_ID,CHR_POS_19)%>%summarise(`P-VALUE`=min(`P-VALUE`)))
	z=merge(x,y,by=c('CHR_ID','CHR_POS_19','P-VALUE'))
	return(z)
}


pick_closest_genes=function(gwas,gene,n=10){
	x=gwas[,CHR_ID]
	y=gene[,chr]
	same_chr=outer(x,y,FUN=`==`)

	x=gwas[,CHR_POS_19]
	y=gene[,tss]
	dist=abs(outer(x,y,FUN=`-`))

	dist[!same_chr]=Inf
	dist_rank=t(apply(dist,1,rank))
	top=(dist_rank<=n)

	container=list(nrow(top))
	for (i in 1:nrow(top)){
		x=esi[which(top[i,]),]
		x=cbind(x,gwas[i,])
		container[[i]]=x
	}
	y=Reduce(rbind,container)
	return(y)
}


pick_top_esi=function(x){
	x[,max_esi:=max(esi),by=c('CHR_ID','CHR_POS_19')]
	top=x[esi==max_esi,]
	return(top)
}

make_boxplot=function(genes,background){
	x=unique(genes[,.(gene_id,esi)])
	y=unique(background[,.(gene_id,esi)])
	z=rbind(data.table(esi=x$esi,gwas='CAD'),data.table(esi=y$esi,gwas='background'))
	pval=wilcox.test(x$esi,y$esi)$p.value
	p=ggplot(z,aes(x=gwas,y=esi))+geom_boxplot()+xlab('GWAS')+ylab('ESI')+annotate(geom='text',x=1.5,y=1.2,label=paste0('Rank-sum test\np < ',formatC(pval,digits=3)))
	return(p)
}

get_tissue_color=function(color_file='shared/tissue_color.txt'){
	tissue_color=fread(color_file)[,c(1,5)]
	colnames(tissue_color)=c('tissue','color')
	tissue_color=tissue_color%>%filter(tissue!='SF')
	color=tissue_color$color
	names(color)=as.character(tissue_color$tissue)
	return(color)
}

make_violin_to_compare_ESI_across_tissue=function(x,title){
	color=get_tissue_color()
	x$tissue=reorder(x$tissue,x$esi,FUN=function(x){mean(x,na.rm=T)})
	p=ggplot(x,aes(tissue,esi,fill=tissue))+geom_violin()+scale_fill_manual(guid='none',values=color)+ylab('ESI')+xlab('Tissue')+coord_flip()+ggtitle(title)
	return(p)
}


bootstrap_difference_in_mean=function(x,y,size=1000){
	stopifnot(colnames(x)==colnames(y))
	diff_mat=matrix(nrow=size,ncol=ncol(x))
	colnames(diff_mat)=colnames(x)
	n=nrow(x)
	m=nrow(y)
	for (i in 1:size){
		xb=x[sample(1:n,n,replace=T),]
		yb=y[sample(1:m,m,replace=T),]
		xmean=apply(xb,2,function(x) {mean(x,na.rm=T)})
		ymean=apply(yb,2,function(x) {mean(x,na.rm=T)})
		diff=xmean-ymean
		diff_mat[i,]=diff
	}
	diff_df=as.data.frame(diff_mat)
	diff_long=melt(diff_df,id.var=NULL,value.name='diff_in_mean',variable.name='tissue')
	return(diff_long)
}


main=function(gwas,gene_location,esi_all,trt){
	message('selecting trait genes...')
	trait=gwas[MAPPED_TRAIT==trt,]
	trait=rmdupvar(trait)
	trait_genes=pick_closest_genes(trait,gene_location)

	message('selecting background genes...')
	background=gwas[MAPPED_TRAIT!=trt,]
	background=rmdupvar(background)
	background_genes=pick_closest_genes(background,gene_location)

	message('selecting top genes...')
	top_genes=pick_top_esi(trait_genes)
	top_background=pick_top_esi(background_genes)

	top_genes_all=esi_all[gene%in%unique(top_genes[,gene_id]),]
	top_background_all=esi_all[gene%in%unique(top_background[,gene_id]),]

	message('performing bootstrap...')
	diff=bootstrap_difference_in_mean(top_genes_all[,2:ncol(top_genes_all),with=F],top_background_all[,2:ncol(top_background_all),with=F],size=1000)
	diff$tissue=reorder(diff$tissue,diff$diff_in_mean,function(x) {mean(x,na.rm=T)})
	diff_summary=diff%>%group_by(tissue)%>%summarise(mean=mean(diff_in_mean),lb=quantile(diff_in_mean,0.025),ub=quantile(diff_in_mean,0.975))
	
	message('making plot')
	p=ggplot(diff_summary,aes(tissue,mean))+geom_point()+geom_errorbar(aes(ymin=lb,ymax=ub))+coord_flip()+geom_hline(yintercept=0)
	return(p)
}

calc_diff_with_matched_snp=function(gwas,matched){
	stopifnot(class(gwas)=='character')
	stopifnot('data.table'%in%class(matched))

	x=as.data.table(str_split_fixed(gwas,':',2))
	setnames(x,c('CHR_ID','CHR_POS_19'))
	x[,CHR_ID:=as.integer(CHR_ID)]
	x[,CHR_POS_19:=as.integer(CHR_POS_19)]
	y=pick_closest_genes(x,esi)
	z=pick_top_esi(y)
	w=esi_all[gene%in%unique(z[,gene_id]),]
	w_long=melt(w,id.var='gene',variable.name='tissue',value.name='esi')
	gwas_mean=w_long%>%group_by(tissue)%>%summarise(mean_esi=mean(esi,na.rm=T))

	diff=data.table()
	for (i in 1:ncol(matched)){
		x=as.data.table(str_split_fixed(unlist(matched[,i,with=F]),':',2))
		setnames(x,c('CHR_ID','CHR_POS_19'))
		x[,CHR_ID:=as.integer(CHR_ID)]
		x[,CHR_POS_19:=as.integer(CHR_POS_19)]
		y=pick_closest_genes(x,esi)
		z=pick_top_esi(y)
		w=esi_all[gene%in%unique(z[,gene_id]),]
		w_long=melt(w,id.var='gene',variable.name='tissue',value.name='esi')
		matched_mean=w_long%>%group_by(tissue)%>%summarise(mean_esi=mean(esi,na.rm=T))
		diff=rbind(diff,data.table(tissue=cad_mean$tissue,diff=cad_mean$mean_esi-matched_mean$mean_esi))
	}
	return(diff)
}

# Read NHGRI GWAS catalog:
# gwas=fread('../data/gwas/nhgri_catalog/gwas_catalog_v1.0.1-associations_e87_r2017-01-23.tsv')
# gwas=reformat(gwas)
# gwas=liftOver(gwas)
# gwasp=prune(gwas)
# fwrite(gwasp,'../data/gwas/nhgri_catalog/gwas_catalog_v1.0.1-associations_e87_r2017-01-23.prune.tsv',sep='\t')
gwasp=fread('../data/gwas/nhgri_catalog/gwas_catalog_v1.0.1-associations_e87_r2017-01-23.prune.tsv')
gwasp=gwasp%>%arrange(CHR_ID,CHR_POS_19)


# Get ESI scores:
esi=fread('../processed_data/160715/esi.quant_norm.hcasmc.txt')
bed=fread('../data/gtex/gencode.v19.genes.v6p.patched_contigs_genetypes.bed')%>%setnames(c('chr','start','end','strand','gene_id','type'))
esi=merge(esi,bed,by='gene_id')
esi[,tss:=ifelse(strand=='+',start,end)]
esi=esi%>%arrange(chr,tss)


# Select variant for CAD or CHD:
cad=gwasp[MAPPED_TRAIT=='coronary heart disease'|MAPPED_TRAIT=='coronary artery disease',]
cad=rmdupvar(cad)

# Output CAD SNPs:
fwrite(cad[,.(CHR_ID,CHR_POS_19)],'../processed_data/160715/SNPsnap/cad_gwas_variants.chr_pos',col.names=F,sep=':')

# Select genes around CAD variants:
cad_genes=pick_closest_genes(cad,esi)


# Select background variants: 
background=gwasp[MAPPED_TRAIT!='coronary heart disease'&MAPPED_TRAIT!='coronary artery disease',]
background=rmdupvar(background)


# Select genes around background variants (takes time):
background_genes=pick_closest_genes(background,esi)


# Make boxplot to compare CAD and background genes:
p1=make_boxplot(cad_genes,background_genes)
save_plot('../figures/160715/hcasmc_specific_gene_and_GWAS.quant_norm.pdf',p1)


# Pick the gene with highest ESI for each loci:
top_cad_genes=pick_top_esi(cad_genes)
top_background_genes=pick_top_esi(background_genes)


# Make boxplot to compare top CAD and background genes: 
p2=make_boxplot(top_cad_genes,top_background_genes)
save_plot('../figures/160715/hcasmc_specific_gene_and_GWAS.top_esi.quant_norm.pdf',p2)


# Read ESI for GTEx tissues and HCASMC: 
esi_all=fread('../processed_data/160715/esi.quant_norm.all_tissues.txt') 


# Select genes around CAD variants:
cad_genes_all_tissue=esi_all[gene%in%unique(cad_genes[,gene_id]),]
cad_genes_long=melt(cad_genes_all_tissue,id.var='gene',variable.name='tissue',value.name='esi')


# Make violin plot on ESI distribution across all tissues for top CAD gene (one per CAD loci):
top_cad_genes_all_tissue=esi_all[gene%in%unique(top_cad_genes[,gene_id]),]
top_cad_genes_long=melt(top_cad_genes_all_tissue,id.var='gene',variable.name='tissue',value.name='esi')
p6=make_violin_to_compare_ESI_across_tissue(top_cad_genes_long,title='top CAD-neighboring genes')
save_plot('../figures/160715/hcasmc_specific_gene_and_GWAS.all_tissue.top_cad_genes.quant_norm.pdf',p6,base_width=8,base_height=8)


# Make violin plot on ESI distribution across all tissues for top non-CAD gene (one per GWAS loci):
top_non_cad_genes_all_tissue=esi_all[gene%in%unique(top_background_genes[,gene_id]),]
top_non_cad_genes_long=melt(top_non_cad_genes_all_tissue,id.var='gene',variable.name='tissue',value.name='esi')
p7=make_violin_to_compare_ESI_across_tissue(top_non_cad_genes_long,title='top non CAD-neighboring genes')
save_plot('../figures/160715/hcasmc_specific_gene_and_GWAS.all_tissue.top_non_cad_genes.quant_norm.pdf',p7,base_width=8,base_height=8)


# Caclulate the difference in ESI between top CAD and top background genes: 
diff_long=bootstrap_difference_in_mean(top_cad_genes_all_tissue[,2:ncol(top_cad_genes_all_tissue),with=F],top_non_cad_genes_all_tissue[,2:ncol(top_non_cad_genes_all_tissue),with=F],size=10000)
diff_long$tissue=reorder(diff_long$tissue,diff_long$diff_in_mean,mean)


# Plot the difference:
diff_summary=diff_long%>%group_by(tissue)%>%summarise(mean=mean(diff_in_mean),lb=quantile(diff_in_mean,0.025),ub=quantile(diff_in_mean,0.975))
p=ggplot(diff_summary,aes(tissue,mean))+geom_point()+geom_errorbar(aes(ymin=lb,ymax=ub))+coord_flip()+geom_hline(yintercept=0)
save_plot('../figures/160715/difference_in_ESI_between_top_CAD_and_background_genes.quant_norm.pdf',p,base_width=8,base_height=8)


# p=main(gwasp,esi,esi_all,'Alzheimers disease')


# Make scatterplot of ESI score (for QC): 
pdf('../figures/160715/esi_scatterplots.quant_norm.pdf')
ggplot(esi_all,aes(HCASMC,`Whole Blood`))+geom_point()+geom_abline(intercept=0,slope=1,color='red')
ggplot(esi_all,aes(HCASMC,`Adipose - Subcutaneous`))+geom_point()+geom_abline(intercept=0,slope=1,color='red')
ggplot(esi_all,aes(HCASMC,`Artery - Coronary`))+geom_point()+geom_abline(intercept=0,slope=1,color='red')
ggplot(esi_all,aes(HCASMC,`Lung`))+geom_point()+geom_abline(intercept=0,slope=1,color='red')
dev.off()


# Read SNPsnap output: 
snpsnap=fread('../processed_data/160715/SNPsnap/SNPsnap_output/matched_snps.txt')


# Calculate difference between top CAD and SNPsnap-matched SNPs:  
diff=calc_diff_with_matched_snp(snpsnap[,Input_SNP],snpsnap[,2:1001,with=F])
diff$tissue=reorder(diff$tissue,diff$diff,median)
diff_summary=res%>%group_by(tissue)%>%summarise(mean=mean(diff),lb=quantile(diff,0.025),ub=quantile(diff,0.975))


# Plot the difference between top CAD and SNPsnap-matched SNPs:
p=ggplot(diff_summary,aes(tissue,mean))+geom_point()+geom_errorbar(aes(ymin=lb,ymax=ub))+coord_flip()+geom_hline(yintercept=0)
save_plot('../figures/160715/difference_in_ESI_between_top_CAD_and_matched_snp.quant_norm.pdf',p,base_width=8,base_height=8)



