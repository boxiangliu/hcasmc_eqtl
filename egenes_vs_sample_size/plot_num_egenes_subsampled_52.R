# library
library(data.table)
library(dtplyr)
library(dplyr)
library(cowplot)
library(stringr)


# read file:
num_egenes=fread('../processed_data/egenes_vs_sample_size/num_egenes_subsampled_52.txt')


# make tissue names pretty: 
num_egenes$tissue=str_replace(num_egenes$tissue,'_',' - ')%>%str_replace_all('_',' ')


# turn tissue into factor, sort by number of egenes 
num_egenes=num_egenes%>%arrange(num_egenes)
num_egenes$tissue=factor(num_egenes$tissue,levels=num_egenes$tissue)


# label artery, heart and smooth muscle: 
num_egenes$Type=NULL
num_egenes$Type='Other'
num_egenes$Type=ifelse(str_detect(num_egenes$tissue,'HCASMC'),'HCASMC',num_egenes$Type)
num_egenes$Type=ifelse(str_detect(num_egenes$tissue,'Artery'),'Artery',num_egenes$Type)
num_egenes$Type=ifelse(str_detect(num_egenes$tissue,'fibroblasts'),'Fibroblast',num_egenes$Type)
num_egenes$Type=ifelse(str_detect(num_egenes$tissue,'Heart'),'Heart',num_egenes$Type)
num_egenes$Type=ifelse(str_detect(num_egenes$tissue,'Colon|Uterus|Esophagus|Stomach|Lung|Vagina|Intestine'),'Smooth Muscle',num_egenes$Type)
num_egenes$Type=factor(num_egenes$Type,levels=c('HCASMC','Artery','Fibroblast','Heart','Smooth Muscle','Other'))


# plot:
num_egenes$size=ifelse(num_egenes$Type=='Other',1,2)
hcasmc_egenes=num_egenes[tissue=='HCASMC',num_egenes]
# p=ggplot(num_egenes,aes(x=num_egenes,y=tissue,color=Type,label=num_egenes))+geom_point()+theme_bw()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+xlab('Number of eGenes')+ylab('Tissue')+scale_color_manual(values=c('red','tan4','linen','deeppink','steelblue','black'))+annotate('text',x=600,y='HCASMC',label=hcasmc_egenes)
p=ggplot(num_egenes,aes(x=num_egenes,y=tissue,color=Type,label=num_egenes,size=size))+geom_point()+theme_bw()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+xlab('Number of eGenes')+ylab('Tissue')+scale_color_manual(values=c(2:6,1))+scale_size(range=c(1,3),guide=F)
save_plot('../figures/egenes_vs_sample_size/num_egenes_subsampled_52.pdf',p,base_height=6)


# get the number of egenes:
range(num_egenes[Type=='Artery',num_egenes])
range(num_egenes[Type=='Smooth Muscle',num_egenes])
range(num_egenes[Type=='Fibroblast',num_egenes])
range(num_egenes[Type=='Heart',num_egenes])

