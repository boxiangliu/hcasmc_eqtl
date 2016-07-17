# command args: 
args=commandArgs(T,T)
ld_file=args$ld
bim_file=args$bim
variants_file=args$variants
figure_file=args$figure

# ld_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/LD2/target_variants.ld'
# bim_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/LD2/target_variants.bim'
# variants_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/compare_hcasmc_and_gtex2/combined.txt'
# figure_file='/srv/persistent/bliu2/HCASMC_eQTL/figures/160615/ld_distribution_by_tissue.pdf'


# read ld file
ld=fread(ld_file,header=F)
bim=fread(bim_file)
variants=fread(variants_file,header=F)


# create id from chr and pos:
bim[,id:=paste(V1,V4,sep='_')]


# cast ld into a matrix:
ld=as.matrix(ld)
colnames(ld)=bim$id
rownames(ld)=bim$id


# create id from chr and pos: 
setnames(variants,c('chr','pos','gene_id','variant_id','type'))
variants[,id:=paste(chr,pos,sep='_')]


# subset LD matrix to GWAS columns: 
gwas=ld[,colnames(ld)%in%variants[type=='GWAS',id]]


# calcualte max LD per row/eQTL: 
max_ld=c()
for (i in 1:nrow(gwas)){
	max_ld=c(max_ld,max(gwas[i,],na.rm=T))
}
gwas_max_ld=data.table(id=rownames(gwas),max_ld=max_ld)


# calculate LD for each tissue:
tissue_max_ld=data.table()
for (t in unique(variants$type)){
	if (t=='GWAS') next()
	idx=gwas_max_ld$id%in%variants[type==t,id]
	temp=data.table(type=t,gwas_max_ld[idx,])
	tissue_max_ld=rbind(tissue_max_ld,temp)
}
tissue_max_ld[,num:=.N,by=type]


# make boxplot for ld:
p=ggplot(tissue_max_ld,aes(x=type,y=max_ld,label=num))+geom_boxplot()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+xlab('Tissue')+ylab('LD')+background_grid(major='x')+geom_text(aes(y=0.5),angle=90,vjust=0.5)
save_plot(figure_file,p,base_height=7, base_width=7)

# blowup: 
p=ggplot(tissue_max_ld,aes(x=type,y=max_ld,label=num))+geom_boxplot()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+xlab('Tissue')+ylab('LD')+background_grid(major='x')+geom_text(aes(y=0.5),angle=90,vjust=0.5)+ylim(0,0.4)
save_plot(str_replace(figure_file,'pdf','blowup.pdf'),p,base_height=7, base_width=7)
