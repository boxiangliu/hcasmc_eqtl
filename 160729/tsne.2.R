library(Rtsne) # Load package
library(RColorBrewer)
library(foreach)
library(doMC)
registerDoMC(cores=30)

# read rpkm: 
rpkm=fread('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160729/combined.rpkm',header=T)


# read column data: 
col_data=fread('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160729/combined.col',header=T)
col_data$tissue=as.factor(col_data$tissue)

# read tissue names: 
tissue_names_file='/srv/persistent/bliu2/HCASMC_eQTL/scripts/160603/collapsed_tissue_names.txt'
tissue_names=fread(tissue_names_file,header=T)


# collapse sub-tissues to reduce the number of tissues:
idx=match(col_data$tissue,tissue_names$original)
col_data$collapsed=tissue_names$collapsed[idx]


# cast rpkm into matrix:
row_names=rpkm$Name
rpkm[,Name:=NULL]
rpkm=as.matrix(rpkm)
rownames(rpkm)=row_names 


# transpose rpkm and normalized the mean of each column (gene) to zero:
rpkm2=scale(t(rpkm),center=T,scale=F)


# make sure the rows are unique: 
rpkm3=unique(rpkm2)


# sanity check; rpkm3 and col_data should be consistent
stopifnot(rownames(rpkm3)==col_data$sample) 


tsne_out_list=foreach(seed = 1:1000) %dopar% {
	# perform t-SNE:
	set.seed(seed)
	# tsne_out=Rtsne(rpkm3)
	tsne_out2=Rtsne(rpkm3,check_duplicates=F,pca=F)


	# extract coordinates from t-SNE output: 
	tsne_Y=data.table(tsne_out2$Y)
	tissue=as.factor(col_data$collapsed)
	tsne_Y$tissue=tissue
	setnames(tsne_Y,c('x','y','tissue'))


	# find centroid of each tissue:
	tsne_Y[,x_centroid:=median(x),by='tissue']
	tsne_Y[,y_centroid:=median(y),by='tissue']
	tsne_Y[,tissue_abb:=toupper(substr(tissue,1,3))]


	# prepare color palette:
	color=as.numeric(tissue)
	qualitative_color=brewer.pal(12, "Paired")
	pal=colorRampPalette(qualitative_color)
	pal=pal(length(levels(tissue)))


	# make scatter plot: 
	pdf(sprintf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160729/tsne.seed%s.pdf',seed),height=8,width=8)
	par(mar=c(8.1,4.1,4.1,7.1))
	par(xpd=T) # turn of clipping to allow legend outside of plotting area
	plot(tsne_Y$x,tsne_Y$y,col=pal[color],pch=19,xlab='Coordinate 1',ylab='Coordinate 2')
	legend(60,60,levels(tissue),col=pal[unique(color)],pch=19,bty='n')
	text(x=tsne_Y$x_centroid,y=tsne_Y$y_centroid,labels=tsne_Y$tissue_abb,cex=0.7)
	dev.off()

	tsne_out2
}


cost=sapply(tsne_out_list, function(x) {y=x$itercost; y[length(y)]})
which.min(cost) # 25
