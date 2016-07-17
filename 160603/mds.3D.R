#!/usr/bin/env Rscript 
# boxiang liu
# durga
# perform multidimentional scaling 
# example: 
# Rscript hclust.R -rpkm=<rpkm_file> -coldata=<coldata_file> -figure=<out_figure>

# library:
library('R.utils')
library('MASS')
library('RColorBrewer')
library('plot3D')
library('rgl')
attach(mtcars)
plot3d(wt, disp, mpg, col="red", size=3)

# command line arguments: 
args=commandArgs(T,T)
rpkm_file=args$rpkm
coldata_file=args$coldata
tissue_names_file=args$tissue_names
figure_dir=args$figure
# rpkm_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined.rpkm'
# coldata_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined.col'
# tissue_names_file='/srv/persistent/bliu2/HCASMC_eQTL/scripts/160603/collapsed_tissue_names.txt'
# figure_dir='/srv/persistent/bliu2/HCASMC_eQTL/figures/160603/'


# read input: 
rpkm=fread(rpkm_file,header=T)
col_data=fread(coldata_file,header=T)
tissue_names=fread(tissue_names_file,header=T)


# create rpkm matrix:
row_data=rpkm$Name
rpkm_mat=as.matrix(subset(rpkm,select=-Name))
rownames(rpkm_mat)=row_data


# calculate distance based on pearson correlation: 
pearson=cor(rpkm_mat)
dist=as.dist(1-pearson)


# perform multidimensional scaling: 
mds_res=isoMDS(dist,k=3)


# collapse sub-tissues to reduce the number of tissues:
idx=match(col_data$tissue,tissue_names$original)
col_data$collapsed=tissue_names$collapsed[idx]


# find centroid of each tissue: 
mds=data.table(mds_res$points,keep.rownames=T)
colnames(mds)=c('sample','x','y','z')
tissue=factor(col_data$collapsed,levels=unique(col_data$collapsed))
mds$tissue=tissue
mds[,x_centroid:=median(x),by='tissue']
mds[,y_centroid:=median(y),by='tissue']
mds[,z_centroid:=median(z),by='tissue']
mds[,tissue_abb:=toupper(substr(tissue,1,3))]


# prepare color palette: 
color=as.numeric(tissue)
qualitative_color=brewer.pal(12, "Paired")
pal=colorRampPalette(qualitative_color)
pal=pal(length(levels(tissue)))


# make scatter plot:
plot3d(mds$x,mds$y,mds$z,col=pal[color],xlab='Coordinate 1',ylab='Coordinate 2',zlab='Coordinate 3')
text3d(mds$x_centroid,mds$y_centroid,mds$z_centroid,texts=mds$tissue_abb)
legend3d('topright',levels(tissue),col=pal[unique(color)],pch=19,bty='n')
# movie3d(spin3d(axis=c(0,1,0)),duration=10,dir='../figures/160603/',convert=F)

