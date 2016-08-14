# Make heatmap based on metasoft m-value: 
library(gplots)


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
