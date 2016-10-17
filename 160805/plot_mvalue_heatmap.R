# Make heatmap based on metasoft m-value: 
library(gplots)
library(data.table)
library(dplyr)
library(cowplot)
library(stringr)


#### function:
getMvalue=function(x,n){x%>%select((17+n):(17+2*n-1))}
getPvalue=function(x,n){x%>%select(17:(17+n-1))}


#### main: 
# read metasoft output: 
metasoft=fread('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/metasoft_output/metasoft_output.1_22.mcmc.txt',skip=1)
metasoft[,V107:=NULL]


# get mvalue: 
n_tissue=45
mvalue=getMvalue(metasoft,n_tissue)


# add column names:
study_name=unlist(fread('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/Metasoft_tissue_order.alphabetical.txt',header=F))
setnames(mvalue,study_name)


# calculate correlation: 
C=cor(mvalue,use='pairwise.complete.obs') # or use pairwise.complete.obs


# # hierarchical clustering: 
# ord=hclust(as.dist(correlation))$order


# # melt correlation matrix C:
# C1=as.data.frame(C)
# C1$tissue1=rownames(C1)
# Cm=melt(C1,id.vars='tissue1',variable.name='tissue2',value.name='correlation')


# # set factor level: 
# Cm$tissue1=factor(Cm$tissue1,levels=rownames(C)[ord],labels=str_replace_all(rownames(C)[ord],"_"," "))
# Cm$tissue2=factor(Cm$tissue2,levels=rownames(C)[ord],labels=str_replace_all(rownames(C)[ord],"_"," "))


# # make heatmap:
# ggplot(Cm,aes(tissue1,tissue2))+geom_tile(aes(fill=correlation))+scale_fill_gradient2(low='blue',high='red')

# create color palette: 
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 100)

# make heatmap: 
pdf('../figures/160805/heatmap_by_mvalue_correlation.pdf',width=10,height=8)
heatmap.2(C,trace='none',col=my_palette,srtCol=45,margins=c(10,15))
dev.off()


# plot HCASMC strip: 
to_plot=data.frame(cor=C[row.names(C)=='HCASMC'],tissue=row.names(C))
to_plot$Type=NULL
to_plot$Type='Other'
to_plot$Type=ifelse(str_detect(to_plot$tissue,'HCASMC'),'HCASMC',to_plot$Type)
to_plot$Type=ifelse(str_detect(to_plot$tissue,'Artery'),'Artery',to_plot$Type)
to_plot$Type=ifelse(str_detect(to_plot$tissue,'fibroblasts'),'Fibroblast',to_plot$Type)
to_plot$Type=ifelse(str_detect(to_plot$tissue,'Heart'),'Heart',to_plot$Type)
to_plot$Type=ifelse(str_detect(to_plot$tissue,'Colon|Uterus|Esophagus|Stomach|Lung|Vagina|Intestine'),'Smooth Muscle',to_plot$Type)
to_plot$Type=factor(to_plot$Type,levels=c('HCASMC','Artery','Fibroblast','Heart','Smooth Muscle','Other'))
to_plot$size=ifelse(to_plot$Type=='Other',1,2)


p=ggplot(to_plot,aes(y=reorder(tissue,cor,FUN=mean),x=cor,color=Type,size=size))+geom_point()+theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1))+xlab('Correlation')+ylab('Tissue')+scale_color_manual(values=c(2:6,1))+scale_size(range=c(1,3),guide=F)
save_plot('../figures/160805/hcasmc_mvalue_correlation.pdf',p,base_height=6)

# perform Wilcox rank sum test: 
to_plot$rank=rank(-to_plot$cor)
artery=to_plot[to_plot$Type=='Artery','rank']
smooth_muscle=to_plot[to_plot$Type=='Smooth Muscle','rank']
heart=to_plot[to_plot$Type=='Heart','rank']
fibroblast=to_plot[to_plot$Type=='Fibroblast','rank']
other=to_plot[to_plot$Type=='Other','rank']

wilcox.test(artery,other_than_artery,alternative='less') # p-value = 0.2794
wilcox.test(heart,other,alternative='less') # p-value = 0.6966
wilcox.test(smooth_muscle,other,alternative='less') # p-value = 0.02369
wilcox.test(fibroblast,other,alternative='less') # p-value = 0.8276
