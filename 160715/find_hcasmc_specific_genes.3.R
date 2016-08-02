# library:
library(DESeq2)
library(gap)
library(gplots)
source('/srv/persistent/bliu2/HCASMC_eQTL/scripts/utils.R')


# read col_data:
load('../processed_data/160715/find_housekeeping_genes.RData')
col_data=as.data.frame(colData(dds)[,c('sample','tissue')])


# construct seq expression set:
residuals=counts(dds,normalized=T)


# get ranks of residuals: 
ranks=getRank(-residuals,dimension=1) # takes a minute


# for each gene, get the lowest rank for HCASMC samples for that gene:
idx=which(col_data$tissue=='HCASMC')
rank_max=apply(ranks[,idx],1,max)


# get p-values:
N=nrow(col_data)
n=length(idx)
pvalue=getPvalue(rank_max,N,n)


# permutation pvalues:
M=length(rank_max)
random_rank_max=numeric(length=M)
set.seed(2)
for (i in 1:M) random_rank_max[i]=max(sample(1:510,10))
perm_pvalue=getPvalue(random_rank_max,N,n)



# make qqplot for p-values:
pdf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/qqplot_pvalues.2.pdf')
qqunif(pvalue)
to_plot=qqunif(perm_pvalue,plot.it=F)
points(to_plot)
legend('topleft',legend=c('nominal','permuted'),col=c('blue','black'),pch=19)
dev.off()


# get FDR:
padjust=p.adjust(pvalue,method='BH')
cutoff=1e-3
sum(padjust<cutoff) 


# subset to HCASMC specific genes with FDR 5%:
hcasmc_genes=residuals[rownames(residuals)%in%names(padjust[padjust<cutoff]),]


# perform hierarchical clustering on HCASMC specific genes: 
pearson=cor(t(hcasmc_genes))
dist=as.dist(1-pearson)
hr <- hclust(dist, method="complete")


# assign clusters: 
# mycl <- cutree(hr, h=max(hr$height/1.2))
mycl <- cutree(hr, h=max(hr$height/1.15))
max(mycl) # 7

# get a color palette equal to the number of clusters
clusterCols <- rainbow(length(unique(mycl)))


# create vector of colors for row side bar:
myClusterSideBar <- clusterCols[mycl]


# create vector of colors for column side bar: 
tissue_colors=ifelse(col_data$tissue[1:100]=='HCASMC',"#FF0000","#0000FF")


# choose a color palette for the heat map
myheatcol <- rev(redgreen(75))


# draw heatmap:
pdf('../figures/160715/heatmap.3.pdf')
heatmap.2(hcasmc_genes[,1:100], main="Hierarchical Cluster", Rowv=as.dendrogram(hr), Colv=T, dendrogram="both", scale="row", col=myheatcol, density.info="none", trace="none", RowSideColors= myClusterSideBar, ColSideColors=tissue_colors)
dev.off()


# read row data: 
row_data=fread('../processed_data/160715/combined.count',header=T)%>%dplyr::select(Name,Description)


# cast cluster assignment into a data frame: 
mycl1=data.frame(gene_id=names(mycl),cluster=mycl)
mycl1=merge(mycl1,row_data,by.x='gene_id',by.y='Name')
mycl1=mycl1%>%dplyr::rename(gene_name=Description)
mycl1=mycl1%>%arrange(cluster)


# output each cluster to a different file: 
for (i in 1:max(mycl)) {
	# subset to genes in cluster i: 
	cluster=mycl1[mycl1$cluster==i,'gene_name']

	# remove number after dot from ENSEMBL ID: 
	cluster=str_split_fixed(cluster,'\\.',n=2)[,1]

	# output cluster 1: 
	write.table(cluster,file=sprintf('../processed_data/160715/cluster%s.txt',i),quote=F,row.names=F,col.names=F)
}


# performed panther analysis using each cluster:
# no code. 
