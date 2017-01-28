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


# Read NHGRI GWAS catalog:
gwas=fread('../data/gwas/nhgri_catalog/gwas_catalog_v1.0.1-associations_e87_r2017-01-23.tsv')
gwas=reformat(gwas)
gwas=liftOver(gwas)
Reading liftover chains
Mapping coordinates
gwasp=prune(gwas)
gwasp=gwasp%>%arrange(CHR_ID,CHR_POS_19)


# Get ESI scores:
esi=fread('../processed_data/160715/esi.hcasmc.txt')
bed=fread('../data/gtex/gencode.v19.genes.v6p.patched_contigs_genetypes.bed')%>%setnames(c('chr','start','end','strand','gene_id','type'))
esi=merge(esi,bed,by='gene_id')
esi[,tss:=ifelse(strand=='+',start,end)]
esi=esi%>%arrange(chr,tss)


# Select variant for CAD or CHD:
cad=gwasp[MAPPED_TRAIT=='coronary heart disease'|MAPPED_TRAIT=='coronary artery disease',]
cad=rmdupvar(cad)


# Select genes around CAD variants:
cad_genes=pick_closest_genes(cad,esi)


# Select variants for all other GWAS: 
background=gwasp[MAPPED_TRAIT!='coronary heart disease'&MAPPED_TRAIT!='coronary artery disease',]
background=rmdupvar(background)


# Select genes around background variants:
background_genes=pick_closest_genes(background,esi)


# Make boxplot to compare CAD and background genes:
p1=make_boxplot(cad_genes,background_genes)
save_plot('../figures/160715/hcasmc_specific_gene_and_GWAS.pdf',p1)


# Pick the gene with highest ESI for each loci:
top_cad_genes=pick_top_esi(cad_genes)
top_background_genes=pick_top_esi(background_genes)


# Make boxplot to compare CAD and background genes: 
p2=make_boxplot(top_cad_genes,top_background_genes)
save_plot('../figures/160715/hcasmc_specific_gene_and_GWAS.top_esi.pdf',p2)


# Read ESI for GTEx tissues and HCASMC: 
esi_all=fread('../processed_data/160715/esi.all_tissues.txt') 


# Select genes around CAD variants:
cad_genes_all_tissue=esi_all[gene%in%unique(cad_genes[,gene_id]),]
cad_genes_long=melt(cad_genes_all_tissue,id.var='gene',variable.name='tissue',value.name='esi')


# Make violin plot to compare ESI distribution across all tissue: 
make_violin_to_compare_ESI_across_tissue=function(){
	tissue_color=fread('shared/tissue_color.txt')[,c(1,5)]
	colnames(tissue_color)=c('tissue','color')
	tissue_color=tissue_color%>%filter(tissue!='SF')
	color=tissue_color$color
	names(color)=as.character(tissue_color$tissue)
	cad_genes_long=merge(cad_genes_long,tissue_color,by='tissue')
	cad_genes_long$tissue=reorder(cad_genes_long$tissue,cad_genes_long$esi,FUN=function(x){mean(x,na.rm=T)})
	ggplot(cad_genes_long,aes(tissue,esi,fill=tissue))+geom_violin()+scale_fill_manual(guid='none',values=color)+ylab('ESI')+xlab('Tissue')+coord_flip()
}
p3=make_violin_to_compare_ESI_across_tissue()
save_plot('../figures/160715/hcasmc_specific_gene_and_GWAS.all_tissue.pdf',p3,base_width=8,base_height=8)