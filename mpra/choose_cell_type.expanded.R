library(data.table)
library(dplyr)
library(dtplyr)
library(stringr)
library(gtools)
library(gplots)
library(cowplot)

# Functions:
filter_metadata=function(x){
	y=x%>%filter(Assembly=='hg19',`Audit ERROR`=='')
	y[,max_size:=max(Size),by=`Biosample term name`]
	z=y[Size==max_size,.(`File accession`,`Biosample term name`,`Biosample type`)]
	z[,`File accession`:=paste0(`File accession`,'.bed.gz')]
	setnames(z,c('file','epigenome','type'))
	return(z)
}

sort_bed_by_chrom_pos=function(x){
	x=as.data.frame(x)
	x$sortcol=paste(x[,1],x[,2],sep='_')
	ord=mixedorder(x[,'sortcol'])
	y=x[ord,]
	y$sortcol=NULL
	y=as.data.table(y)
}

calculate_jaccard=function(annotation,dir){
	jaccard=matrix(NA,nrow=length(annotation$file),ncol=length(annotation$file))
	colnames(jaccard)=rownames(jaccard)=annotation$epigenome
	for (i in 1:nrow(jaccard)){
		for (j in i:ncol(jaccard)){
			file1=paste(dir,str_replace(annotation$file[i],'.gz',''),sep='/')
			file2=paste(dir,str_replace(annotation$file[j],'.gz',''),sep='/')
			res=system(sprintf('bedtools jaccard -g ../../shared/genomes/hg19.chrom.sorted.sizes -a %s -b %s',file1,file2),intern = T)
			jaccard[i,j]=as.numeric(str_split_fixed(res[2],'\\t',4)[3])
		}
	}
	jaccard[lower.tri(jaccard)]=t(jaccard)[lower.tri(jaccard)]
	return(jaccard)
}

select_tissue=function(jaccard,tissue){
	x=jaccard[rownames(jaccard)==tissue,]
	x=data.frame(epigenome=names(x),jaccard=x)%>%filter(epigenome!=tissue)
	x$epigenome=reorder(x$epigenome,x$jaccard,mean)
	return(x)
}

gwas2bed=function(gwas,window_size){
	gwas[,start:=pos-window_size-1]
	gwas[,end:=pos+window_size]
	if (!str_detect(gwas[1,chr],'chr')) {gwas[,chr:=paste0('chr',chr)]}
	gwas=gwas[,.(chr,start,end,markername)]
	return(gwas)
}

extract_basename=function(x){
	y=basename(x)
	y=str_replace(y,'.bed.gz','')
	y=str_split_fixed(y,'-|_',2)[,1]
	return(y)
}



# Create temporary directory: 
tmp_dir=tempdir()
if (!dir.exists(tmp_dir)) {dir.create(tmp_dir)}

# Rename sample by cell type: 
annotation=fread('../data/encode/dnase_seq/metadata.tsv')%>%filter(Assembly=='hg19')
annotation=filter_metadata(annotation)
annotation=rbind(data.frame(file='2305_ppr.IDR0.1.filt.narrowPeak.gz',epigenome='HCASMC',type='primary cell'),annotation)


# Process DHS data:
out_dir=paste(tmp_dir,'top50000',sep='/')
if (!dir.exists(out_dir)) dir.create(out_dir)
for (i in annotation$file){
	x=fread(paste('zcat ../processed_data/mpra/DHS_expanded/',i,sep='/'))
	setnames(x,c('chrom','chromStart','chromEnd','name','score','strand','signalValue','pValue','qValue','peak'))
	y=x[chrom%in%paste0('chr',1:22),] # keeping only autosomes.
	y=y%>%arrange(desc(signalValue)) # sort peaks by signalValue
	if (nrow(y)>=50000){
		z=y[1:50000,] # take the top 50,000 peaks 
	} else {
		message(sprintf('%s has less than 50,000 peaks',i))
		z=y
	}
	z=sort_bed_by_chrom_pos(z)
	fwrite(z,paste(out_dir,str_replace(i,'.gz',''),sep='/'),sep='\t',col.names=F)
}


# Calculate the Jaccard index:
jaccard=calculate_jaccard(annotation,out_dir)


# Make heatmap of jaccard index: 
pdf('../figures/mpra/epigenome_jaccard_index_heatmap.expanded.pdf',width=20,height=20)
heatmap.2(jaccard,trace='none',margin=c(10,10))
dev.off()


# Plot Jaccard index against HCASMC:
jaccard_hcasmc=select_tissue(jaccard,'HCASMC')
pdf('../figures/mpra/epigenome_jaccard_index_against_hcasmc.expanded.pdf',width=8,height=24)
ggplot(jaccard_hcasmc,aes(epigenome,jaccard))+geom_point()+coord_flip()+ylab('Jaccard index')+xlab('')+background_grid(major = "xy", minor = "y")
dev.off()


# Read in GWAS loci:
gwas=fread('../processed_data/mpra/rAggr/maf0.01_dist500kb_rsq0.8_1kgp3.uniq.txt')


# Expand each GWAS variant up- and downstream 500 bases, then convert to bed format: 
gwas=gwas2bed(gwas,500)


# Output GWAS regions: 
gwas_outfile=paste(tmp_dir,'gwas_regions.bed',sep='/')
fwrite(gwas,gwas_outfile,col.names=F,sep='\t')


# Overlap GWAS regions with DHS regions:
in_dir=paste(tmp_dir,'top50000',sep='/')
out_dir=paste(tmp_dir,'top50000_gwas',sep='/')
if (!dir.exists(out_dir)) dir.create(out_dir)
for (i in annotation$file){
	in_file=paste(in_dir,str_replace(i,'.gz',''),sep='/')
	out_file=paste(out_dir,str_replace(i,'.gz',''),sep='/')
	system(sprintf('bedtools intersect -wa -a %s -b %s > %s',in_file,gwas_outfile,out_file))
}

# Calculate jaccard index: 
jaccard2=calculate_jaccard(annotation,out_dir)


# Make heatmap of jaccard index: 
pdf('../figures/mpra/epigenome_jaccard_index_heatmap.gwas_regions.expanded.pdf',width=20,height=20)
heatmap.2(jaccard2,trace='none',margin=c(10,10))
dev.off()

# Plot Jaccard index against HCASMC:
jaccard2_hcasmc=select_tissue(jaccard2,'HCASMC')
pdf('../figures/mpra/epigenome_jaccard_index_against_hcasmc.gwas_regions.expanded.pdf',width=8,height=24)
ggplot(jaccard2_hcasmc,aes(epigenome,jaccard))+geom_point()+coord_flip()+ylab('Jaccard index')+xlab('')+background_grid(major = "xy", minor = "y")
dev.off()


# Remove temporary directory: 
unlink(tmp_dir,recursive=T)