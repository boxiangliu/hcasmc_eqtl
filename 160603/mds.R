#!/usr/bin/env Rscript 
# boxiang liu
# durga
# perform multidimentional scaling 
# example: 
# Rscript hclust.R -rpkm=<rpkm_file> -coldata=<coldata_file> -figure=<out_figure>

# library:
library(dplyr)
library(data.table)
library('R.utils')
library('MASS')
library('RColorBrewer')
library(cowplot)
library(ggrepel)

# command line arguments: 
args=commandArgs(T,T)
rpkm_file=args$rpkm
coldata_file=args$coldata
tissue_names_file=args$tissue_names
figure_path=args$figure
# rpkm_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined.rpkm'
# coldata_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined.col'
# tissue_names_file='/srv/persistent/bliu2/HCASMC_eQTL/scripts/160603/collapsed_tissue_names.2.txt'
# figure_path='/srv/persistent/bliu2/HCASMC_eQTL/figures/160603/'


# read input: 
rpkm=fread(rpkm_file,header=T)
col_data=fread(coldata_file,header=T)
tissue_names=fread(tissue_names_file,header=T)


# create rpkm matrix:
row_data=rpkm$Name
rpkm_mat=as.matrix(subset(rpkm,select=-Name))
rownames(rpkm_mat)=row_data


# calculate distance based on pearson correlation:
message('calculating distance...')
pearson=cor(rpkm_mat)
dist=as.dist(1-pearson)


# perform multidimensional scaling: 
message('multidimensional scaling...')
mds_res=isoMDS(dist,k=2)


# collapse sub-tissues to reduce the number of tissues:
idx=match(col_data$tissue,tissue_names$original)
col_data$collapsed=tissue_names$collapsed[idx]


# find centroid of each tissue: 
mds=data.table(mds_res$points,keep.rownames=T)
colnames(mds)=c('sample','x','y')
tissue=factor(col_data$collapsed,levels=unique(col_data$collapsed))
mds$tissue=tissue
mds[,x_centroid:=median(x),by='tissue']
mds[,y_centroid:=median(y),by='tissue']
mds[,tissue_abb:=toupper(substr(tissue,1,3))]


# prepare color palette: 
color=as.numeric(tissue)
qualitative_color=brewer.pal(12, "Paired")
pal=colorRampPalette(qualitative_color)
pal=pal(length(levels(tissue)))


# make scatter plot:
message('plotting...')
set.seed(42)
mds_centroid=unique(mds[,.(tissue,x_centroid,y_centroid)])
p1=ggplot(mds,aes(x=x,y=y,color=tissue))+geom_point()+geom_text_repel(data=mds_centroid,aes(x=x_centroid,y=y_centroid,label=tissue),color='black')+theme_bw()+xlab('MDS Coordinate 1')+ylab('MDS Coordinate 2')+theme(axis.title=element_text(size=15))
save_plot(paste0(figure_path,'/mds.pdf'),p1,base_height=6,base_width=8.5)
p2=ggplot(mds,aes(x=x,y=y,color=tissue))+geom_point()+xlim(-0.1,0.1)+ylim(-0,0.2)+geom_text_repel(data=mds_centroid,aes(x=x_centroid,y=y_centroid,label=tissue),size=7,force=3,color='black')+theme_bw()+scale_color_discrete(guide=F)+theme(axis.text=element_blank(),axis.title=element_blank())
save_plot(paste0(figure_path,'/mds.cropped.pdf'),p2)
save.image('../processed_data/160603/mds.Rdata')

