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


# get q-values:
padjust=p.adjust(pvalue,method='BH')
sum(padjust<0.05) # 6684


# output pvalues:
temp=read.table('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/combined.count',header=T,check.names=F)%>%dplyr::select(gene_id=Name,gene_name=Description)
output=data.frame(gene_id=names(pvalue),pvalue=pvalue,padjust=padjust)
output=merge(output,temp,by='gene_id')
setcolorder(output,c('gene_id','gene_name','pvalue','padjust'))
output=output%>%arrange(pvalue)



# synthetic SMC specific genes: 
output[str_detect(output$gene_id,'ENSG00000122786'),] # CALD1
output[str_detect(output$gene_id,'ENSG00000026025'),] # VIM
output[str_detect(output$gene_id,'ENSG00000133026'),] # MYH10
output[str_detect(output$gene_id,'ENSG00000167460'),] # TPM4
output[str_detect(output$gene_id,'ENSG00000114115'),] # RBP1
# > output[str_detect(output$gene_id,'ENSG00000122786'),] # CALD1
#                 gene_id gene_name       pvalue      padjust
# 3523 ENSG00000122786.15     CALD1 9.766996e-06 0.0001316586
# > output[str_detect(output$gene_id,'ENSG00000026025'),] # VIM
#                gene_id gene_name       pvalue      padjust
# 1982 ENSG00000026025.9       VIM 7.126712e-09 1.691097e-07
# > output[str_detect(output$gene_id,'ENSG00000133026'),] # MYH10
#                gene_id gene_name       pvalue     padjust
# 2141 ENSG00000133026.8     MYH10 1.503665e-08 3.33892e-07
# > output[str_detect(output$gene_id,'ENSG00000167460'),] # TPM4
#                gene_id gene_name       pvalue      padjust
# 480 ENSG00000167460.10      TPM4 2.154211e-15 2.055935e-13
# > output[str_detect(output$gene_id,'ENSG00000114115'),] # RBP1
#                 gene_id gene_name    pvalue padjust
# 26405 ENSG00000114115.5      RBP1 0.6277038       1
write.table(output,file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/hcasmc_specific_genes.2.txt',quote=F,row.names=F,sep='\t')


# output the rank for GSEA:
set.seed(2)
ranks=data.frame(gene_name=output$gene_name,rank=rank(-output$padjust,ties.method='random'))
write.table(ranks,file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/hcasmc_specific_genes.rank.2.rnk',quote=F,row.names=F,sep='\t')
write.table(ranks%>%dplyr::select(gene_name),file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/hcasmc_specific_genes.rank.2.gene_name_only.rnk',quote=F,row.names=F,sep='\t')

# make heatmap of gene expression profiles: 
# head(output)
# to_plot=residuals[rownames(residuals)%in%output$gene_id[1:200],1:50]
# head(to_plot)
# output[2,]
# plot(unlist(residuals[rownames(residuals)=='ENSG00000150261.2',]))
# output[1,]
# rownames(residuals)=='ENSG00000123500.5'
# heatmap.2(as.matrix(to_plot),Rowv=F,Colv=F,dendrogram='none',trace='none')



# select genes highly expressed in HCASMC: 
selected=output[output$padjust<1e-15,'gene_id']
residuals_selected=residuals[rownames(residuals)%in%selected,]
residuals_selected=residuals_selected[,1:100]
tissue_colors=ifelse(col_data$tissue[1:100]=='HCASMC',"#FF0000","#0000FF")
pdf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/heatmap.2.pdf')
heatmap.2(residuals_selected, col=topo.colors(75), scale="row", ColSideColors=tissue_colors,
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
dev.off()


