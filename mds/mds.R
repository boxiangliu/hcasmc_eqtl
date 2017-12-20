#!/usr/bin/env Rscript 
# boxiang liu
# durga
# perform multidimentional scaling 
# example: 
# Rscript hclust.R -rpkm=<rpkm_file> -coldata=<coldata_file> -figure=<out_figure>

# library:
library(dplyr)
library(dtplyr)
library(data.table)
library(R.utils)
library(MASS)
library(RColorBrewer)
library(cowplot)
library(ggrepel)

get_color_map=function(){
	color=fread('shared/tissue_color.txt')
	color[,tissue_color_hex:=max(tissue_color_hex),by=tissue]
	color_map=color$tissue_color_hex
	names(color_map)=color$tissue
	return(color_map)
}


color_file='shared/tissue_color.txt'
rpkm_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/mds/combine_hcasmc_and_gtex_rpkm/combined.rpkm'
coldata_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/mds/combine_hcasmc_and_gtex_rpkm/combined.col'
rowdata_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/mds/combine_hcasmc_and_gtex_rpkm/combined.row'
figure_path='/srv/persistent/bliu2/HCASMC_eQTL/figures/mds/'
out_path='../processed_data/mds/mds/'
if (!dir.exists(figure_path)){dir.create(figure_path,recursive=TRUE)}
if (!dir.exists(out_path)){dir.create(out_path,recursive=TRUE)}

# read input: 
rpkm=fread(rpkm_file,header=T)
col_data=fread(coldata_file,header=T)
row_data=fread(rowdata_file,header=T)


# create rpkm matrix:
rpkm_mat=as.matrix(subset(rpkm,select=-Name))
rownames(rpkm_mat)=row_data$Name


# calculate distance based on pearson correlation:
message('calculating distance...')
pearson=cor(rpkm_mat)
dist=as.dist(1-pearson)


# perform multidimensional scaling: 
message('multidimensional scaling...')
mds_res=isoMDS(dist,k=2)


# collapse sub-tissues to reduce the number of tissues:
tissue_color=fread('shared/tissue_color.txt')[,list(tissue_site_detail,tissue_group=tissue)]
col_data=merge(col_data,tissue_color,by.x='tissue',by.y='tissue_site_detail')


# find centroid of each tissue: 
mds=data.table(mds_res$points,keep.rownames=T)
colnames(mds)=c('sample_id','x','y')
tissue=factor(col_data$tissue_group,levels=unique(col_data$tissue_group))
mds=merge(mds,col_data,by='sample_id')
mds[,x_centroid:=median(x),by='tissue_group']
mds[,y_centroid:=median(y),by='tissue_group']


# prepare color palette:
color_map=get_color_map()


# make scatter plot:
message('plotting...')
set.seed(42)
mds_centroid=unique(mds[,list(tissue_group,x_centroid,y_centroid)])
mds_centroid[,label:=ifelse(tissue_group %in% c('Liver','Muscle',
	'HCASMC','Fibroblasts','Heart','Brain','Blood','LCL',
	'Testis','Pancreas','Spleen','Pituitary'),as.character(tissue_group),'')]

xmin=-0.05
xmax=0.075
ymin=0
ymax=0.1

p1=ggplot(mds,aes(x=x,y=y,color=tissue_group))+
	geom_point()+scale_color_manual(values=color_map)+
	geom_text(data=mds_centroid,aes(x=x_centroid,
		y=y_centroid,label=label),color='black',
		fontface=ifelse(mds_centroid$label=='HCASMC','bold','plain'))+
	theme_bw()+xlab('MDS Coordinate 1')+ylab('MDS Coordinate 2')+
	theme(axis.title=element_text(size=15),
		legend.position = c(0.18,0.7),
		legend.title=element_blank(),
		legend.box.margin=margin(0,0,0,0),
		legend.background=element_rect(fill=alpha('white',0)))+
	geom_rect(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,color='red',linetype=2,alpha=0)


save_plot(paste0(figure_path,'/mds.pdf'),p1,base_height=8,base_width=8.5)

p2=ggplot(mds,aes(x=x,y=y,color=tissue_group))+
	geom_point(alpha=ifelse(mds$tissue_group%in%c(
		'Adipose','Breast','Bladder',
		'Blood Vessel','Cervix','Colon','Esophagus','Fibroblasts',
		'Fallopian Tube','Heart','Vagina','Uterus'),1,0))+
	xlim(xmin,xmax)+
	ylim(ymin,ymax)+
	geom_text_repel(data=mds_centroid,
		aes(x=x_centroid,y=y_centroid,label=tissue_group),
		fontface=ifelse(mds_centroid$label=='HCASMC','bold','plain'),
		force=3,color='black')+
	theme_bw()+
	scale_color_manual(guide='none',values=color_map)+
	theme(axis.text=element_blank(),
		axis.title=element_blank(),
		axis.ticks=element_blank())

save_plot(paste0(figure_path,'/mds.cropped.pdf'),p2)

p3=ggdraw()+draw_plot(p1)+draw_plot(p2,0.65,0.5,0.33,0.48)
save_plot(sprintf('%s/mds.inset.pdf',figure_path),p3, base_aspect_ratio = 1.0, base_height=8)

save.image(sprintf('%s/mds.Rdata',out_path))

