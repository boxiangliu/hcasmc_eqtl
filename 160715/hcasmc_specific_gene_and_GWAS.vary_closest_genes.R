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

rm_overlap_var=function(x,y,min_dist=3e6){
	### remove variants in x with less than min_dist from any variants in y
	stopifnot(class(x[,CHR_ID])==class(y[,CHR_ID]))
	stopifnot(class(x[,CHR_POS_19])==class(y[,CHR_POS_19]))
	same_chr=outer(x[,CHR_ID],y[,CHR_ID],`==`)
	dist=abs(outer(x[,CHR_POS_19],y[,CHR_POS_19],`-`))
	dist[!same_chr]=Inf
	dist_to_closest=apply(dist,1,min)
	z=x[dist_to_closest>=min_dist]
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

plot_distance=function(x,title){
	lb=quantile(x[,dist],probs=c(0.005))
	ub=quantile(x[,dist],probs=c(0.995))
	hist(x$dist,breaks=50,xlab='Distance',main=title)
	abline(v=lb,col='red',lty=2)
	abline(v=ub,col='red',lty=2)
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


main=function(gwas,gene_location,esi_all,trait_variants,background_variants,n_genes,title){
	message('selecting trait genes...')
	trait_genes=pick_closest_genes(trait_variants,gene_location,n=n_genes)

	message('selecting background genes...')
	background_genes=pick_closest_genes(background_variants,gene_location,n=n_genes)

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
	color=get_tissue_color() # get tissue color
	p=ggplot(diff_summary,aes(tissue,mean))+geom_errorbar(aes(ymin=lb,ymax=ub,color=tissue),width=0.1,size=1)+geom_point(size=1)+coord_flip()+geom_hline(yintercept=0)+ggtitle(title)+scale_color_manual(guide='none',values=color)
	return(p)
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


# Select genes around CAD variants:
cad_genes=pick_closest_genes(cad,esi,n=10)


# Select background variants: 
background=gwasp[MAPPED_TRAIT!='coronary heart disease'&MAPPED_TRAIT!='coronary artery disease',]
background=rmdupvar(background)


# Remove variants from background set that overlap the CAD set:
background=rm_overlap_var(background,cad)


# Select genes around background variants (takes time):
background_genes=pick_closest_genes(background,esi,n=10)


# Pick the gene with highest ESI for each loci:
top_cad_genes=pick_top_esi(cad_genes)
top_background_genes=pick_top_esi(background_genes)


# Read ESI for GTEx tissues and HCASMC: 
esi_all=fread('../processed_data/160715/esi.quant_norm.all_tissues.txt') 


# Plot the difference in ESI between top CAD and top background genes using various numbers of closest genes: 
pdf('../figures/160715/difference_in_ESI_between_top_CAD_and_background_genes.vary_closest_genes.pdf')
for (n in c(3,5,10,15,20)){
	print(main(gwasp,esi,esi_all,cad,background,n_genes=n,paste(n,'closest gene')))
}
dev.off()



