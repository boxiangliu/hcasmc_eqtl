library(data.table)
library(dplyr)
library(dtplyr)
library(MASS)
library(stringr)
library(cowplot)
library(gplots)
library(doMC)
library(foreach)
registerDoMC(40)

bed_dir='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_filt/'
sorted_bed_dir='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_filt/sorted/'
out_dir='../processed_data/hcasmc_specific_open_chromatin/mds_open_chromatin_merged_peaks/'
fig_dir='../figures/hcasmc_specific_open_chromatin/mds_open_chromatin_merged_peaks/'
if (!dir.exists(sorted_bed_dir)) dir.create(sorted_bed_dir)
if (!dir.exists(out_dir)) dir.create(out_dir)
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

sort_bed=function(annotation,in_dir,out_dir){
	for (i in 1:length(annotation$file)){
		print(sprintf('INFO - %s',annotation$file[i]))
		in_fn=sprintf('%s/%s',in_dir,annotation$file[i])
		out_fn=sprintf('%s/%s',out_dir,annotation$file[i])
		cmd=sprintf('sort -k1,1 -k2,2g "%s" > "%s"',in_fn,out_fn)
		system(cmd)
	}
}

calculate_jaccard=function(annotation,dir){
	jaccard=matrix(NA,nrow=length(annotation$file),ncol=length(annotation$file))
	colnames(jaccard)=rownames(jaccard)=annotation$epigenome
	for (i in 1:nrow(jaccard)){
		# for(j in i:ncol(jaccard)){
		# 	file1=paste(dir,str_replace(annotation$file[i],'.gz',''),sep='/')
		# 	file2=paste(dir,str_replace(annotation$file[j],'.gz',''),sep='/')
		# 	print(sprintf("%s and %s",file1,file2))
		# 	res=system(sprintf('bedtools jaccard -a "%s" -b "%s"',file1,file2),intern = T)
		# 	jaccard[i,j]=as.numeric(str_split_fixed(res[2],'\\t',4)[3])
		# }
		
		jaccard[i,i:ncol(jaccard)]=foreach(j = i:ncol(jaccard),.combine='c')%dopar%{
			file1=paste(dir,str_replace(annotation$file[i],'.gz',''),sep='/')
			file2=paste(dir,str_replace(annotation$file[j],'.gz',''),sep='/')
			print(sprintf("%s and %s",file1,file2))
			res=system(sprintf('bedtools jaccard -a "%s" -b "%s"',file1,file2),intern = T)
			as.numeric(str_split_fixed(res[2],'\\t',4)[3])
		}
	}
	jaccard[lower.tri(jaccard)]=t(jaccard)[lower.tri(jaccard)]
	return(jaccard)
}


detect_duplicates=function(x){
	duplicates=data.table(which(x==1,arr.ind=TRUE),keep.rownames=TRUE)
	duplicates=duplicates[row!=col,]
	pairs_done=c()
	remove=c()
	for (i in seq(nrow(duplicates))){
		pair=unlist(duplicates[i,list(row,col)])
		item=pair['row']
		if (!all(pair%in%pairs_done)){
			remove[length(remove)+1]=item
			pairs_done[(length(pairs_done)+1):(length(pairs_done)+2)]=pair
		}
	}
	return(remove)
}

#----------- Using Jaccard index to measure distance ------------#
annotation=data.table(file=list.files(bed_dir,pattern='.bed'))
annotation[,epigenome:=str_replace(file,'.bed','')]
annotation
sort_bed(annotation,bed_dir,sorted_bed_dir)

jaccard=calculate_jaccard(annotation,sorted_bed_dir)
saveRDS(jaccard,sprintf('%s/jaccard.rds',out_dir))


# Caclulate distance matrix for hierarchical clustering:
dist=as.dist(1-jaccard)




# Multi-dimensional scaling:
mds=isoMDS(dist,k=2)
to_plot=data.table(mds$points,keep.rownames=TRUE)
to_plot[,tissue:=str_split_fixed(rn,'\\.',2)[,1]]
p=ggplot(to_plot,aes(x=V1,y=V2,color=tissue))+geom_point()
save_plot(sprintf('%s/mds_jaccard.pdf',fig_dir),p,base_width=25)


# Plot hierarchical clustering: 
hc=hclust(dist)
pdf(sprintf('%s/hierarchical_clustering_jaccard.pdf',fig_dir),height=50,width=50)
plot(hc)
dev.off()


# Heatmap: 
pdf(sprintf('%s/heatmap.pdf',fig_dir),height=50,width=50)
heatmap.2(jaccard,trace="none",Rowv = FALSE, Colv=FALSE, dendrogram='none')
heatmap.2(jaccard,trace="none")
dev.off()

