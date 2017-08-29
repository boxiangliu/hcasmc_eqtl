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
color_file='shared/tissue_color.txt'

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
mds$tissue=col_data$tissue
mds$collapsed=col_data$collapsed
mds[,x_centroid:=median(x),by='collapsed']
mds[,y_centroid:=median(y),by='collapsed']



# prepare color palette: 
color=fread(color_file,select=c(3,5),col.names=c('tissue','color'))
color=merge(color,tissue_names,by.x='tissue',by.y='original')
color_map=color$color
names(color_map)=color$collapsed


# make scatter plot:
message('plotting...')
set.seed(42)
mds_centroid=unique(mds[,list(collapsed,x_centroid,y_centroid)])
mds_centroid[,label:=ifelse(collapsed %in% c('Liver','Muscle',
	'HCASMC','Fibroblasts','Heart','Brain','Blood','LCL',
	'Testis','Pancreas','Spleen','Pituitary'),as.character(collapsed),'')]

p1=ggplot(mds,aes(x=x,y=y,color=collapsed))+
	geom_point()+scale_color_manual(values=color_map)+
	geom_text(data=mds_centroid,aes(x=x_centroid,
		y=y_centroid,label=label),color='black',
		fontface=ifelse(mds_centroid$label=='HCASMC','bold','plain'))+
	theme_bw()+xlab('MDS Coordinate 1')+ylab('MDS Coordinate 2')+
	theme(axis.title=element_text(size=15),
		legend.position = c(0.18,0.73),
		legend.title=element_blank(),
		legend.box.margin=margin(0,0,0,0),
		legend.background=element_rect(fill=alpha('white',0)))


save_plot(paste0(figure_path,'/mds.pdf'),p1,base_height=6,base_width=8.5)

p2=ggplot(mds,aes(x=x,y=y,color=collapsed))+
	geom_point(alpha=ifelse(mds$collapsed%in%c('HCASMC','Muscle',
		'Fibroblasts','Skin','Artery','Heart',
		'Esophagus','Adipose','Vagina','Colon',
		'Uterus','Nerve'),1,0))+xlim(-0.1,0.1)+
	ylim(-0,0.2)+
	geom_text_repel(data=mds_centroid,
		aes(x=x_centroid,y=y_centroid,label=collapsed),
		fontface=ifelse(mds_centroid$label=='HCASMC','bold','plain'),
		force=3,color='black')+
	theme_bw()+
	scale_color_manual(guide='none',values=color_map)+
	theme(axis.text=element_blank(),
		axis.title=element_blank())

save_plot(paste0(figure_path,'/mds.cropped.pdf'),p2)

p3=ggdraw()+draw_plot(p1+theme(axis.title.x=element_text(hjust=0.75)),0,0,1,1)+draw_plot(p2,0,0,0.55,0.4)
save_plot(sprintf('%s/mds.inset.pdf',figure_path),p3, base_aspect_ratio = 1.0, base_height=7)

save.image('../processed_data/160603/mds.Rdata')
# load('../processed_data/160603/mds.Rdata')
