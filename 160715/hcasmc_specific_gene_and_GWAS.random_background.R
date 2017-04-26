library(data.table)
library(dtplyr)
library(dplyr)
library(stringr)
library(cowplot)


# Variables: 
fig_dir='../figures/160715/hcasmc_specific_gene_and_GWAS.random_background/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir)}

# Function:
reformat=function(x){
	if (class(x$CHR_POS)!='numeric') {x[,CHR_POS:=as.integer(CHR_POS)]}
	y=x[!is.na(CHR_POS),]
	y=y[CHR_ID!='X',]
	return(y)
}

liftOver=function(x,chain='/srv/persistent/bliu2/shared/chain_files/hg38ToHg19.over.chain.gz'){
	x[,rowid:=1:nrow(x)]
	x[,CHR_POS:=as.integer(CHR_POS)]
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


# Read NHGRI GWAS catalog:
gwas=fread('../data/gwas/nhgri_catalog/gwas_catalog_v1.0.1-associations_e87_r2017-01-23.tsv')


# Select variants from Nikpay et al:
cad=gwas[`FIRST AUTHOR`=='Nikpay M']
cad=liftOver(cad)
cad=rmdupvar(cad)
cad=cad[,.(CHR_ID,CHR_POS_19,`P-VALUE`,PUBMEDID,`FIRST AUTHOR`,`DISEASE/TRAIT`)]
cad[,snpID:=paste(CHR_ID,CHR_POS_19,sep=':')]


# Read SNPsnap database:
snpsnap=fread('/srv/persistent/bliu2/shared/SNPsnap/ld0.5_collection.tab',select=c(1,2,3,4,5,6,25))
snpsnap[,chr:=str_split_fixed(snpID,':',2)[,1]]
snpsnap=snpsnap[chr%in%c(1:22),]


# Select background variants:
background_ls=list()
for (i in 1:nrow(cad)){

	snpid=unlist(cad[i,snpID])
	print(sprintf('INFO - %s',snpid))
	tmp=snpsnap[snpID==snpid,list(freq_bin,dist_nearest_gene_snpsnap,friends_ld05)]
	freq=tmp$freq_bin
	dist=tmp$dist_nearest_gene_snpsnap
	ld05=tmp$friends_ld05

	step=0
	tmp=data.frame()
	while (nrow(tmp)<500){
		step=step+1
		ub=1+0.1*step
		lb=1-0.1*step
		tmp=snpsnap[(freq_bin<=freq+step)&(freq_bin>=freq-step)&(dist_nearest_gene_snpsnap<=ub*dist)&(dist_nearest_gene_snpsnap>=lb*dist)&(friends_ld05<=ub*ld05)&(friends_ld05>=lb*ld05),]
	}
	print(sprintf('INFO - tolerance: %s',step))
	set.seed(42)
	tmp=tmp[snpID!=snpid]
	background=tmp[sample(1:nrow(tmp),500),]
	background_ls[[snpid]]=background
}

# Read gene annotation: 
annotation=fread('../data/gtex/gencode.v19.genes.v6p.patched_contigs_genetypes.bed',col.names=c('chr','start','end','strand','gene_id','type'))
annotation=annotation[type%in%c('lincRNA','protein_coding')&chr%in%c(1:22),]


# Read expression specificity index:
esi=fread('../processed_data/160715/esi.hcasmc.txt')
esi=merge(esi,annotation,by='gene_id')


# Read expression ranks:
expression_rank=read.table('../processed_data/160715/expression_median_and_rank/rank.tsv',header=T,row.names=1,sep='\t',check.names=F)
tissue='HCASMC'


# For each variants, pick genes in 1Mb window: 
var=data.table(snpID=names(background_ls[1]),status='gwas')
var=rbind(var,data.table(snpID=background_ls[[1]][,snpID],status='background'))
var[,chr:=str_split_fixed(snpID,":",2)[,1]]
var[,pos:=as.integer(str_split_fixed(snpID,":",2)[,2])]

pdf(sprintf('%s/max_rank.pdf',fig_dir))
window=c(1e4,1e5,1e6)
for (w in window){
	var[,c('start','end'):=list(pos-w,pos+w)]
	setkey(var,chr,start,end)
	setkey(annotation,chr,start,end)
	overlap=foverlaps(var,annotation,nomatch=0)

	container=list()
	for (v in unique(overlap$snpID)){
		print(sprintf('INFO - %s',v))
		gene_id=overlap[snpID==v,gene_id]
		status=unique(overlap[snpID==v,status])
		tmp=max(expression_rank[which(rownames(expression_rank)%in%gene_id),tissue])
		container[[v]]=data.frame(rank=tmp,snpID=v,status)
	}

	x=Reduce(rbind,container)
	print(ggplot(x,aes(snpID,rank,color=status))+geom_point()+ggtitle(w))
}
dev.off()

