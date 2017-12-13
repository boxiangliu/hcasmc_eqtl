# Make locuscatter plots
# Boxiang Liu
# 2017-12-07
library(R.utils)
library(data.table)
library(cowplot)
library(ggrepel)
library(stringr)

args=commandArgs(T,T,
	defaults=list(marker_col1='rsid',
				  pval_col1='pval',
				  marker_col2='rsid',
				  pval_col2='pval',
				  title1='Study 1',
				  title2='Study 2',
				  snp=NULL))

in_fn1=args$f1
marker_col1=args$marker_col1
pval_col1=args$pval_col1
title1=args$title1

in_fn2=args$f2
marker_col2=args$marker_col2
pval_col2=args$pval_col2
title2=args$title2

tmp_dir=args$tmp_dir
snp=args$snp

in_fn1='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/rasqual/output_pval/chr15/ENSG00000182511.7_FES.pval.txt'
marker_col1='rsid'
pval_col1='pval'

in_fn2='/srv/persistent/bliu2/HCASMC_eQTL/data//gwas/ukbb/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz'
marker_col2='snptestid'
pval_col2='p-value_gc'

out_dir='../processed_data/finemap/locuszoom/locuscatter/'
if (!dir.exists(out_dir)){dir.create(out_dir,recursive=TRUE)}

read_metal=function(in_fn,marker_col,pval_col){
	if (grepl('.gz',in_fn)){
		d=fread(sprintf('gunzip -c %s',in_fn))
	} else {
		d=fread(in_fn)
	}
	setnames(d,c(marker_col,pval_col),c('rsid','pval'))
	d[,list(rsid,pval,logp=-log10(pval))]
}

get_chr=function(eqtl_fn){
	as.integer(str_replace(unique(fread(eqtl_fn)$chr),'chr',''))
}

extract_population=function(population,
	out_file,
	panel='/srv/persistent/bliu2/shared/1000genomes/phase3v5a/integrated_call_samples_v3.20130502.ALL.panel'){
	panel=fread(panel)
	x=panel[super_pop==population,list(sample,sample)]
	fwrite(x,out_file,sep='\t',col.names=FALSE)
}


subset_vcf=function(vcf_in,rsid,population,vcf_out){
	pop_fn=sprintf('%s/%s.txt',out_dir,pop)
	rsid_fn=sprintf('%s/rsid.txt',out_dir)

	extract_population('EUR',pop_fn)
	write.table(rsid,rsid_fn,sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

	command=sprintf('plink --vcf %s --keep-allele-order --keep %s --extract %s --recode vcf-iid --out %s',vcf_in,pop_fn,rsid_fn,vcf_out)
	print(command)
	system(command)
}


calc_LD=function(rsid,chr,pop,out_dir,
	vcf_dir='/srv/persistent/bliu2/shared/1000genomes/phase3v5a/'){
	pop_fn=sprintf('%s/%s.txt',out_dir,pop)
	rsid_fn=sprintf('%s/rsid.txt',out_dir)

	extract_population('EUR',pop_fn)
	write.table(rsid,rsid_fn,sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

	subset_vcf_prefix=sprintf('%s/%s',out_dir,pop)
	subset_vcf_fn=sprintf('%s/%s.vcf',out_dir,pop)
	subset_vcf(sprintf('%s/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',vcf_dir,chr),
		rsid,population,subset_vcf_prefix)

	command=sprintf('plink --vcf %s --keep-allele-order --r2 --out %s/%s',subset_vcf_fn,out_dir,pop)
	print(command)
	system(command)

	ld=fread(sprintf('%s/%s.ld',out_dir,pop))
	return(ld)
}

assign_color=function(rsid,snp,ld){
	all_snps=unique(ld$SNP_A)

	color_dt=ld[SNP_A==snp,list(rsid=SNP_B,color=cut(R2,breaks=c(0,0.2,0.4,0.6,0.8,1),
			labels=c('navyblue','skyblue','darkgreen','orange','red'),
			include.lowest=TRUE))]
	color_dt=rbind(color_dt,data.table(rsid=all_snps[!all_snps%in%color_dt$rsid],color='blue4'))
	color_dt=rbind(color_dt,data.table(rsid=rsid[!rsid%in%all_snps],color='grey'))
	color_dt[rsid==snp,color:='purple']
	color=as.character(color_dt$color)
	names(color)=color_dt$rsid
	return(color)
}

make_combined_plot=function(merged,title1,title2,ld,snp=NULL){
	if (is.null(snp)){
		snp=merged[which.min(pval1+pval2),rsid]
	} else {
		if(!snp%in%merged$rsid){
			stop(sprintf('%s not found in %s',snp,in_fn1))
		}
	}
	print(sprintf('INFO - %s',snp))

	color=assign_color(merged$rsid,snp,ld)
	shape=ifelse(merged$rsid==snp,23,21)
	names(shape)=merged$rsid
	size=ifelse(merged$rsid==snp,3,2)
	names(size)=merged$rsid
	merged[,label:=ifelse(rsid==snp,rsid,'')]

	p1=make_locuscatter(merged,title1,title2,ld,color,shape,size)
	p2=make_locuszoom(merged[,list(rsid,logp=logp1,label)],title1,ld,color,shape,size)
	p2=p2+theme(axis.text.x=element_blank(),axis.title.x=element_blank())
	p3=make_locuszoom(merged[,list(rsid,logp=logp2,label)],title2,ld,color,shape,size)
	p4=plot_grid(p2,p3,align='v',nrow=2)
	p5=plot_grid(p1,p4)
	return(p5)
}

make_locuscatter=function(merged,title1,title2,ld,color,shape,size){
	p=ggplot(merged,aes(logp1,logp2))+
		geom_point(aes(fill=rsid,size=rsid,shape=rsid),alpha=0.8)+
		geom_point(data=merged[label!=''],aes(logp1,logp2,fill=rsid,size=rsid,shape=rsid))+
		xlab(title1)+ylab(title2)+
		scale_fill_manual(values=color,guide='none')+
		scale_shape_manual(values=shape,guide='none')+
		scale_size_manual(values=size,guide='none')+
		geom_text_repel(aes(label=label))
	return(p)
}
make_locuszoom=function(metal,title,ld,color,shape,size){
	data=merge(metal,unique(ld[,list(chr=CHR_A,pos=BP_A,rsid=SNP_A)]),by='rsid')
	chr=unique(data$chr)
	ggplot(data,aes(x=pos,logp))+
		geom_point(aes(fill=rsid,size=rsid,shape=rsid),alpha=0.8)+
		geom_point(data=data[label!=''],aes(x=pos,logp,fill=rsid,size=rsid,shape=rsid))+
		scale_fill_manual(values=color,guide='none')+
		scale_shape_manual(values=shape,guide='none')+
		scale_size_manual(values=size,guide='none')+
		geom_text_repel(aes(label=label))+
		xlab(paste0('chr',chr))+
		ylab(paste(title,'\n-log10(P)'))
}


d1=read_metal(in_fn1,marker_col1,pval_col1)
d2=read_metal(in_fn2,marker_col2,pval_col2)
chr=get_chr(in_fn1)
merged=merge(d1,d2,by='rsid',suffixes=c('1','2'),all=FALSE)
ld=calc_LD(merged$rsid,chr,'EUR',out_dir)
p=make_combined_plot(merged,'eQTL','UKBB',ld,'rs2521501')
pdf('1.pdf',height=4,width=8)
print(p)
dev.off()