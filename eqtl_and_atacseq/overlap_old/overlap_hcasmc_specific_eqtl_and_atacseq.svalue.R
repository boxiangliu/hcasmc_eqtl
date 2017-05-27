# Calculate eQTL and ATACseq overlap. eQTLs are stratified by s-value
# defined as p(HCASMC)*\prod{(1-p(GTEx tissue))}
# where p is the m-value returned by metasoft. 

## Library: 
library(data.table)
library(dplyr)
library(stringr)
library(cowplot)
library(ggrepel)
library(dtplyr)
library(Hmisc)

## Function: 
read_metasoft=function(in_file,tissue_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/Metasoft_tissue_order.alphabetical.txt'){
	header=scan(in_file,what='character',nlines=1,sep='\t')
	metasoft=fread(paste0("awk '{if (NF==106) {print $0}}' ", in_file),skip=1)
	metasoft[,V107:=NULL]


	# read tissue name (study name):
	study_name=unlist(fread(tissue_file,header=F))


	# set metasoft column names:
	col_names=c(paste('pvalue',study_name,sep='_'),paste('mvalue',study_name,sep='_'))
	col_names=c(header[1:16],col_names)
	stopifnot(ncol(metasoft)==length(col_names))
	setnames(metasoft,col_names)
	return(metasoft)
}


get_mvalue=function(x){
	mvalue=x%>%dplyr::select(contains('mvalue'))%>%as.data.frame()
	rownames(mvalue)=x$RSID
	colnames(mvalue)=str_replace(colnames(mvalue),'mvalue_','')
	return(mvalue)
}


calc_svalue=function(x,toi){
	x[x==0]=1e-6
	x[x==1]=1-1e-6
	y=data.frame(log(x[,toi]))
	z=log(1-x[,colnames(x)!=toi])
	w=cbind(y,z)
	s=data.frame(s=apply(w,1,sum))
	rownames(s)=rownames(x)
	return(s)
}


define_peak_region=function(x,size=100){
	x[,peak_pos:=start+peak]
	x[,start:=as.integer(peak_pos-size)]
	x[,end:=as.integer(peak_pos+size)]
	x[,peak_pos:=NULL]
	return(x)
}

subset2bestQTL=function(x,by,rank){
	setnames(x,c(by,rank),c('fid','logpval'))
	x=x%>%group_by(fid)%>%mutate(is_best=(logpval==max(logpval)))
	x=x%>%filter(is_best==TRUE)
	setnames(x,c('fid','logpval'),c(by,rank))
	x[,is_best:=NULL]
	return(as.data.table(x))
}

get_gene_name_and_id=function(gencode_file){
	gencode=fread(gencode_file)
	gencode=gencode%>%filter(V3=="gene")
	gene_id=str_extract(gencode$V9,'(?<=gene_id ")(ENSG.+?)(?=";)')
	gene_name=str_extract(gencode$V9,'(?<=gene_name ")(.+?)(?=";)')
	stopifnot(length(gene_id)==length(gene_name))
	x=data.table(gene_id=gene_id,gene_name=gene_name)
	return(x)
}

name2id=function(name,gene_name_and_id){
	gene_id=c()
	for (n in name){
		x=gene_name_and_id%>%filter(gene_name==n)%>%select(gene_id)%>%unlist()
		if (length(x)>1) {stop("Name is not unique!")}
		gene_id=c(gene_id,x)
	}
	names(gene_id)=NULL
	return(gene_id)
}

id2name=function(id,gene_name_and_id){
	gene_name=gene_name_and_id$gene_name[match(id,gene_name_and_id$gene_id)]
	return(gene_name)
}

parse_svalue_rownames=function(svalue){
	rownames=str_split_fixed(rownames(svalue),'_',6)[,1:5]%>%as.data.table()
	setnames(rownames,c('fid','chr','pos','ref','alt'))
	if (!str_detect(rownames$chr[1],'chr')){
		rownames=rownames%>%mutate(chr=paste0('chr',chr))
	}
	svalue=data.table(cbind(rownames,svalue))
	return(svalue)
}

append_svalue=function(x,svalue,id=c('fid','chr','pos')){
	setnames(x,id,c('fid','chr','pos'))
	if (!is.character(x$pos)) {
		x$pos=as.character(x$pos)
	}
	if (!is.character(svalue$pos)){
		svalue$pos=as.character(svalue$pos)
	}
	x$tmp_id=apply(x%>%select(fid,chr,pos),1,function(x){paste(x,collapse='_')})
	svalue$tmp_id=apply(svalue%>%select(fid,chr,pos),1,function(x){paste(x,collapse='_')})
	x$svalue=svalue[match(x$tmp_id,svalue$tmp_id),s]
	x[,tmp_id:=NULL]
	svalue[,tmp_id:=NULL]
	x$pos=as.integer(x$pos)
	svalue$pos=as.integer(svalue$pos)
	return(x)
}

read_metasoft_dir=function(dir,pattern,svalue_tissue='HCASMC'){
	# Find all files that matches pattern: 
	in_files=list.files(path=dir,pattern=pattern,full.names=T)

	# Initialize list container: 
	container=list()

	for (in_file in in_files){
		# Read metasoft result:
		message(in_file)
		if (str_detect(in_file,'metasoft_output.X')){next()}

		metasoft=read_metasoft(in_file)

		# Subset svalue to bestQTL:
		mvalue=get_mvalue(metasoft)


		# Subset to rows where HCASMC is not NA and fill NAs in other tissues:
		mvalue=mvalue[!is.na(mvalue$HCASMC),]
		mvalue[is.na(mvalue)]=NA_fill


		# Calculate S-value:
		svalue=calc_svalue(mvalue,svalue_tissue)

		# Add to list: 
		container[[length(container)+1]]=svalue
	}

	return(container)
}

plot_overlap=function(overlapped,ratio=1,g=20,plot.it=T,breaks=c(0,1,2,3,4,5,6,Inf)){
	# overlapped[,bin:=cut2(logpval,g=g)]
	overlapped[,bin:=cut(logpval,breaks=breaks)]
	overlapped[,pct:=mean(!is.na(start))/ratio,by=c('bin','specific')]
	pct=overlapped%>%select(bin,specific,pct)%>%unique()
	pct$lb=as.numeric(str_extract(pct$bin,'(?<=[\\(\\[])(.+?)(?=,)'))
	pct=pct%>%arrange(lb)
	if (plot.it){
		p=ggplot(pct,aes(x=lb,y=100*pct))+geom_point(size=5,alpha=0.8)+stat_smooth()+xlab('-log10(P-value)')+ylab('Percentage overlap')
		return(p)
	} else {
		return(pct)
	}
}

## Variables: 
N_TISSUE=45
NA_fill=0.5
atacseq_file='../data/atacseq/fbs/2305/out/peak/idr/optimal_set/2305_ppr.IDR0.1.filt.narrowPeak'
fig_dir='../figures/eqtl_and_atacseq/'
eqtl_file='../data/eQTL/rasqual/expressedGenes.padj.txt'
gencode_file='/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf'


# Read HCASMC eQTL:
eqtl=fread(paste0("cat ", eqtl_file, ' | cut -f1,2,3,4,5,6,7,12,25,26,29'))
eqtl=eqtl[r2_rsnp>0.8,]
eqtl[,logpval:=-log10(pval)]
bestQTL=subset2bestQTL(x=eqtl,by='fid',rank='logpval')
bestQTL=bestQTL%>%mutate(start=pos,end=pos)


# Read metasoft output and calculate svalue:
container=read_metasoft_dir('../processed_data/160805/metasoft_output_subsample_52_p1e-2/','metasoft_output.*.mcmc.txt')

# Concatenate: 
svalue=Reduce(rbind,container)

# Parse rownames of s-value to be fid, chr, pos, ref, alt columns:
svalue=parse_svalue_rownames(svalue)


# Convert the gene ID to gene symbol (name): 
gene_name_and_id=get_gene_name_and_id(gencode_file)
svalue$fid=id2name(svalue$fid,gene_name_and_id)


# Append the s-value to bestQTL dataframe:
bestQTL=append_svalue(bestQTL,svalue)
bestQTL=bestQTL[!is.na(bestQTL$svalue)]


# Calculate s-value quantiles: 
bestQTL$specific=bestQTL[,cut2(svalue,g=5)]


# Read ATACseq data: 
atac_hcasmc=fread(atacseq_file)%>%setnames(.,c('chr','start','end','name','score','strand','signalValue','pValue','qValue','peak'))
atac_hcasmc=define_peak_region(atac_hcasmc)


# Overlap eQTL with ATACseq: 
setkey(bestQTL,chr,start,end)
setkey(atac_hcasmc,chr,start,end)
overlapped=foverlaps(bestQTL,atac_hcasmc)


# Plot eQTL and ATACseq overlap percentage:
pct=plot_overlap(overlapped,plot.it=F)
p1=ggplot(pct,aes(x=lb,y=100*pct,color=specific))+geom_point(size=5,alpha=0.8)+geom_line()+xlab('-log10(P-value)')+ylab('Percentage overlap')+ylim(0,5)
save_plot(paste0('../figures/eqtl_and_atacseq/','eqtl_overlap_open_chromatin_by_specificity.pdf'),p1,base_height=6,base_width=6)
