library(data.table)
library(dplyr)
library(dtplyr)
library(stringr)
library(gtools)
library(gplots)
library(cowplot)

# Functions 
sort_bed_by_chrom_pos=function(x){
	x=as.data.frame(x)
	x$sortcol=paste(x[,1],x[,2],sep='_')
	ord=mixedorder(x[,'sortcol'])
	y=x[ord,]
	y$sortcol=NULL
	y=as.data.table(y)
}

gwas2bed=function(gwas,window_size){
	gwas[,start:=pos-window_size-1]
	gwas[,end:=pos+window_size]
	if (!str_detect(gwas[1,chr],'chr')) {gwas[,chr:=paste0('chr',chr)]}
	gwas=gwas[,.(chr,start,end,markername)]
	return(gwas)
}

# Read DHS file names: 
files=list.files('../processed_data/mpra/DHS','*narrowPeak.gz',full.names=T)
tmp_dir=tempdir()
if (!dir.exists(tmp_dir)) {dir.create(tmp_dir)}


# Process DHS data:
out_files=c()
for (i in files){
	x=fread(paste('zcat', i))
	setnames(x,c('chrom','chromStart','chromEnd','name','score','strand','signalValue','pValue','qValue','peak'))
	y=x[chrom%in%paste0('chr',1:22),] # keeping only autosomes.
	y=y%>%arrange(desc(pValue)) # sort peaks by -log10(pvalue)
	if (nrow(y)>=50000){
		z=y[1:50000,] # take the top 50,000 peaks 
	} else {
		message(sprintf('%s has less than 50,000 peaks',i))
		z=y
	}
	out_file=paste(tmp_dir,str_replace(basename(i),'.gz',''),sep='/')
	out_files=c(out_files,out_file)
	z=sort_bed_by_chrom_pos(z)
	fwrite(z,out_file,sep='\t',col.names=F)
}


# Rename sample by cell type: 
annotation=fread('../data/roadmap/roadmap_sample_annotation.tsv')%>%filter(MARK=='DNase')%>%select(EID,`Standardised epigenome name`)%>%setnames(c('sample','epigenome'))
annotation=rbind(data.frame(sample='2305',epigenome='HCASMC'),annotation)
epigenome=annotation[match(str_split_fixed(basename(out_files),'[_-]',2)[,1],annotation$sample),epigenome]


# Calculate the Jaccard index:
jaccard=matrix(NA,nrow=length(out_files),ncol=length(out_files))
colnames(jaccard)=rownames(jaccard)=epigenome
for (i in 1:length(out_files)){
	for (j in 1:length(out_files)){
		file1=out_files[i]
		file2=out_files[j]
		res=system(sprintf('bedtools jaccard -a %s -b %s',file1,file2),intern = T)
		jaccard[i,j]=as.numeric(str_split_fixed(res[2],'\\t',4)[3])
	}
}


# Make heatmap of jaccard index: 
pdf('../figures/mpra/epigenome_jaccard_index_heatmap.pdf',width=12,height=12)
heatmap.2(jaccard,trace='none',margin=c(20,20))
dev.off()


# Plot Jaccard index against HCASMC:
jaccard_hcasmc=jaccard[rownames(jaccard)=='HCASMC',]
jaccard_hcasmc=data.frame(epigenome=names(jaccard_hcasmc),jaccard=jaccard_hcasmc)%>%filter(epigenome!='HCASMC')
jaccard_hcasmc$epigenome=reorder(jaccard_hcasmc$epigenome,jaccard_hcasmc$jaccard,mean)
pdf('../figures/mpra/epigenome_jaccard_index_against_hcasmc.pdf',width=8,height=8)
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
out_files2=c()
for (i in out_files){
	out_file2=str_replace(i,'narrowPeak','gwasIntersect.narrowPeak')
	out_files2=c(out_files2,out_file2)
	system(sprintf('bedtools intersect -wa -a %s -b %s > %s',i,gwas_outfile,out_file2))
}

jaccard2=matrix(NA,nrow=length(out_files2),ncol=length(out_files2))
colnames(jaccard2)=rownames(jaccard2)=epigenome
for (i in 1:length(out_files2)){
	for (j in 1:length(out_files2)){
		file1=out_files2[i]
		file2=out_files2[j]
		res=system(sprintf('bedtools jaccard -g ../../shared/genomes/hg19.chrom.sorted.sizes -a %s -b %s',file1,file2),intern = T)
		jaccard2[i,j]=as.numeric(str_split_fixed(res[2],'\\t',4)[3])
	}
}

# Make heatmap of jaccard index: 
pdf('../figures/mpra/epigenome_jaccard_index_heatmap.gwas_regions.pdf',width=12,height=12)
heatmap.2(jaccard2,trace='none',margin=c(20,20))
dev.off()

# Plot Jaccard index against HCASMC:
jaccard2_hcasmc=jaccard2[rownames(jaccard2)=='HCASMC',]
jaccard2_hcasmc=data.frame(epigenome=names(jaccard2_hcasmc),jaccard=jaccard2_hcasmc)%>%filter(epigenome!='HCASMC')
jaccard2_hcasmc$epigenome=reorder(jaccard2_hcasmc$epigenome,jaccard2_hcasmc$jaccard,mean)
pdf('../figures/mpra/epigenome_jaccard_index_against_hcasmc.gwas_regions.pdf',width=8,height=8)
ggplot(jaccard2_hcasmc,aes(epigenome,jaccard))+geom_point()+coord_flip()+ylab('Jaccard index')+xlab('')+background_grid(major = "xy", minor = "y")
dev.off()


# Remove temporary directory: 
unlink(tmp_dir,recursive=T)