library(data.table)
library(dplyr)
library(dtplyr)
library(MASS)
library(stringr)

bed_dir='../processed_data/hcasmc_specific_open_chromatin/encode_plus_hcasmc_filt/'
fig_dir='../figures/hcasmc_specific_open_chromatin/mds/'
if (!dir.exists(fig_dir)) dir.create(fig_dir)

# Functions: 
calc_dist=function(sharing){
	dist_mat=matrix(nrow=ncol(sharing),ncol=ncol(sharing))
	colnames(dist_mat)=rownames(dist_mat)=copy(colnames(sharing))
	for (i in 1:ncol(sharing)){
		s1=colnames(sharing)[i]
		setnames(sharing,s1,'sample1')
		for (j in 1:ncol(sharing)){
			if (i==j){
				dist_mat[i,j]=0
			} else {
				s2=colnames(sharing)[j]
				setnames(sharing,s2,'sample2')
				dist_mat[i,j]=sum(sharing[,sample1!=sample2])
				setnames(sharing,'sample2',s2)
			}
		}
		setnames(sharing,'sample1',s1)
	}
	dist=as.dist(dist_mat)
	return(dist)
}


calculate_jaccard=function(annotation,dir){
	jaccard=matrix(NA,nrow=length(annotation$file),ncol=length(annotation$file))
	colnames(jaccard)=rownames(jaccard)=annotation$epigenome
	for (i in 1:nrow(jaccard)){
		for (j in i:ncol(jaccard)){
			file1=paste(dir,str_replace(annotation$file[i],'.gz',''),sep='/')
			file2=paste(dir,str_replace(annotation$file[j],'.gz',''),sep='/')
			print(sprintf("%s and %s",file1,file2))
			res=system(sprintf('bedtools jaccard -a %s -b %s',file1,file2),intern = T)
			jaccard[i,j]=as.numeric(str_split_fixed(res[2],'\\t',4)[3])
		}
	}
	jaccard[lower.tri(jaccard)]=t(jaccard)[lower.tri(jaccard)]
	return(jaccard)
}


#-------- Using peak sharing (indicator variable) to measure distance ---------# 
in_fn='../processed_data/hcasmc_specific_open_chromatin/intersect/intersect.bed'
sharing=fread(in_fn)
sharing[,c('chrom','start','end','num','list'):=list(NULL,NULL,NULL,NULL,NULL)]


# Caclulate distance matrix for hierarchical clustering:
dist=calc_dist(sharing)


# Plot hierarchical clustering: 
hc=hclust(dist)
pdf(sprintf('%s/hierarchical_clustering_peak_sharing.pdf',fig_dir),height=16,width=16)
plot(hc)
dev.off()

# performing PCA:
sharing_t=t(sharing)
pc<-prcomp(sharing_t,center=F,scale.=F,tol=0)


# Plot PCA:
pdf(sprintf('%s/pca_peak_sharing.pdf',fig_dir),height=12,width=12)
plot(pc$x[,1],pc$x[,2])
text(pc$x[,1],pc$x[,2],labels=rownames(pc$x))
screeplot(pc)
dev.off()


# Multi-dimensional scaling:
mds=isoMDS(dist,k=2)
pdf(sprintf('%s/mds_peak_sharing.pdf',fig_dir),width=16,height=16)
with(mds,plot(points[,1],points[,2]))
with(mds,text(points[,1],points[,2],labels=rownames(points)))
dev.off()



#------- Using signal values to measure distance ---------#
in_fn='../processed_data/hcasmc_specific_open_chromatin/intersect/intersect.quant_norm.count'
sharing=fread(in_fn)
setDF(sharing)
rownames(sharing)=sharing$id
sharing$id=NULL


# Caclulate distance matrix for hierarchical clustering:
dist=as.dist(1-cor(sharing,method='pearson'))


# Plot hierarchical clustering: 
hc=hclust(dist)
pdf(sprintf('%s/hierarchical_clustering_signal_value.pdf',fig_dir),height=16,width=16)
plot(hc)
dev.off()



# performing PCA:
sharing_t=t(sharing)
pc<-prcomp(sharing_t,center=T,scale.=T,tol=0)


# Plot PCA:
pdf(sprintf('%s/pca_signal_value.pdf',fig_dir),height=12,width=12)
plot(pc$x[,1],pc$x[,2])
text(pc$x[,1],pc$x[,2],labels=rownames(pc$x))
screeplot(pc)
dev.off()


# Multi-dimensional scaling:
mds=isoMDS(dist,k=2)
pdf(sprintf('%s/mds_signal_value.pdf',fig_dir),width=16,height=16)
with(mds,plot(points[,1],points[,2]))
with(mds,text(points[,1],points[,2],labels=rownames(points)))
dev.off()


#----------- Using Jaccard index to measure distance ------------#
annotation=data.table(file=list.files(bed_dir,pattern='.bed'))
annotation[,epigenome:=str_replace(file,'.bed','')]
jaccard=calculate_jaccard(annotation,bed_dir)


# Caclulate distance matrix for hierarchical clustering:
dist=as.dist(1-jaccard)


# Plot hierarchical clustering: 
hc=hclust(dist)
pdf(sprintf('%s/hierarchical_clustering_jaccard.pdf',fig_dir),height=16,width=16)
plot(hc)
dev.off()


# Multi-dimensional scaling:
mds=isoMDS(dist,k=2)
pdf(sprintf('%s/mds_jaccard.pdf',fig_dir),width=16,height=16)
with(mds,plot(points[,1],points[,2]))
with(mds,text(points[,1],points[,2],labels=rownames(points)))
dev.off()
