# read ld file
ld_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/LD/target_variants.ld'
ld=fread(ld_file,header=F)
bim_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/LD/target_variants.bim'
bim=fread(bim_file)
bim[,id:=paste(V1,V4,sep='_')]

# ld[,id:=paste(CHR_A,BP_A,sep='_')]
# ld[,id2:=paste(CHR_B,BP_B,sep='_')]
ld=as.matrix(ld)
colnames(ld)=bim$id
rownames(ld)=bim$id


variants_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/compare_hcasmc_and_gtex/combined.txt'
variants=fread(variants_file,header=F)

setnames(variants,c('chr','pos','gene_id','variant_id','type'))
variants[,id:=paste(chr,pos,sep='_')]

gwas=ld[,colnames(ld)%in%variants[type=='GWAS',id]]

gwas=ld[colnames(ld)%in%variants[type=='GWAS',id] | ld$id%in%variants[type=='GWAS',id2],]
idx=c()
for (i in 1:nrow(gwas)){
	max_ld=max(gwas[i,],na.rm=T)
	idx=c(idx,max_ld>=0.1)
}
gwas_ld=gwas[idx,]


# calculate the number of variants in LD with GWAS variants:
num_ld_variants=data.table()
for (t in unique(variants$type)){
	if (t=='GWAS') next()
	temp=data.table(type=t,num=sum(variants[type==t,id]%in%rownames(gwas_ld)))
	num_ld_variants=rbind(num_ld_variants,temp)
}


# calcualte the total nubmer of variants:
num_variants=data.table()
for (t in unique(variants$type)){
	if (t=='GWAS') next()
	temp=data.table(type=t,num=sum(variants[type==t,id]%in%rownames(gwas)))
	num_variants=rbind(num_variants,temp)
}



# make plot:
ggplot(num_ld_variants,aes(type,num))+geom_point()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+xlab('Tissue')+ylab('Proportion')+background_grid(major='x')


# calcualte max ld per row: 
max_ld=c()
for (i in 1:nrow(gwas)){
	max_ld=c(max_ld,max(gwas[i,],na.rm=T))
}
gwas_max_ld=data.table(id=rownames(gwas),max_ld=max_ld)


# calculate ld for each tissue:
tissue_max_ld=data.table()
for (t in unique(variants$type)){
	if (t=='GWAS') next()
	idx=gwas_max_ld$id%in%variants[type==t,id]
	temp=data.table(type=t,gwas_max_ld[idx,])
	tissue_max_ld=rbind(tissue_max_ld,temp)
}
tissue_max_ld[,num:=.N,by=type]


# make boxplot for ld:
ggplot(tissue_max_ld,aes(x=type,y=max_ld,label=num))+geom_boxplot()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+xlab('Tissue')+ylab('LD')+background_grid(major='x')+geom_text(aes(y=0.5),angle=90,vjust=0.5)+ylim(0,0.4)


# calculate eQTL LDs:
eqtl_ld=data.table()
for (t in unique(variants$type)){
	if (t=='GWAS') next()
	idx=rownames(ld)%in%variants[type==t,id]
	temp=ld[idx,idx]
	temp=temp[upper.tri(temp)]
	temp=sum(temp>0.8)
	temp=data.table(type=t,num=temp)
	eqtl_ld=rbind(eqtl_ld,temp)
}

ggplot(eqtl_ld,aes(x=type,y=num))+geom_point()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+xlab('Tissue')+ylab('Number of eQTL with LD > 0.8')+background_grid(major='x')


