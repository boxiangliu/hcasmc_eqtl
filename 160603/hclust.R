#!/usr/bin/env Rscript
# boxiang liu
# durga 
# perform hierarchical clustering
# example: 
# Rscript hclust.R -rpkm=<rpkm_file> -coldata=<coldata_file> -figure=<out_figure>

# library:
library('R.utils')


# command line arguments: 
args=commandArgs(T,T)
rpkm_file=args$rpkm
coldata_file=args$coldata
figure_file=args$figure
rpkm_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined.rpkm'
coldata_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined.col'
figure_file='/srv/persistent/bliu2/HCASMC_eQTL/figures/160603/hclust.pdf'

# read input: 
rpkm=fread(rpkm_file,header=T)
col_data=fread(coldata_file,header=T)


# create rpkm matrix:
row_data=rpkm$Name
rpkm_mat=as.matrix(subset(rpkm,select=-Name))
rownames(rpkm_mat)=row_data


# perform clustering: 
pearson=cor(rpkm_mat)
dist=as.dist(1-pearson)
clusters=hclust(dist,method='average')
clusters$labels=ifelse(runif(nrow(col_data))<=0.1,col_data$tissue,'')


# make dendrogram: 
pdf(figure_file,width=200,height=10)
plot(clusters,xlab='Samples', main='Pearson correlation, Average Linkage')
dev.off()


# plot strip colored by tissue:
idx=clusters$order
tissue=col_data$tissue[idx]
tissue=as.factor(tissue)
pal=rainbow(length(levels(tissue)))
pdf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160603/hclust.color.pdf',width=200,height=2)
barplot(rep(1,length(tissue)),col=pal[as.numeric(tissue)],border=NA)
dev.off()


# save environment: 
save.image('../processed_data/160603/hclust.RData')
