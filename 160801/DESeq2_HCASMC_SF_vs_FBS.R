#########
# SETUP
#########
library(dplyr)
library(data.table)
library(DESeq2)
library("BiocParallel")
register(MulticoreParam(10))

#######
# MAIN
#######

####################
# read count matrix
####################
count_wide = fread('../processed_data/160801/rnaseq_dase.combined.count')


# coerce to data.frame: 
count_matrix = count_wide %>% dplyr::select(-Name,-Description) %>% data.frame(check.names=F)
rownames(count_matrix) = count_wide %>% select(Name) %>% unlist()


########################
# create colData matrix
########################
# individual:
individual = str_split_fixed(colnames(count_matrix),'_',2)[,1]

# condition: 
condition = str_split_fixed(colnames(count_matrix),'_',2)[,2] %>% str_replace(.,"_[123]",'')

# batch:
batch2 = c('9090701_sf','8072501_sf','8100901_sf','1051601_sf','9071501_sf','9052004_fbs','9052004_sf')
batch=ifelse(colnames(count_matrix) %in% batch2, 'batch2','batch1')
batch=ifelse(str_detect(colnames(count_matrix),'2305'),'batch3',batch)

# colData: 
colData = data.frame(
  row.names = colnames(count_matrix), 
  condition = condition,
  individual = individual,
  batch = batch
)


########################
# create DESeq2 dataset
########################
dds = DESeqDataSetFromMatrix(countData = count_matrix,
                             colData = colData,
                             design = ~ condition + individual)

# change base level to sf. 
dds$condition = relevel(dds$condition, 'sf') 

# change base level to batch1. 
# dds$batch = relevel(dds$batch, 'batch1')


##############
# run DESeq2 
##############
dds = DESeq(dds, parallel = TRUE)
res = results(dds, contrast = c('condition','fbs','sf'))
resOrdered = res[order(res$padj),]
# sum(resOrdered$padj<1e-5,na.rm=T)


# make MA plot: 
plotMA(res,main='HCASMC FBS vs SF')


# rank genes in descending order based on pvalues: 
ranks = rank(-resOrdered$pvalue,ties.method='random',na.last=F)
temp=data.frame(gene_id=rownames(resOrdered),rank=ranks)
output=merge(temp,count_wide[,.(Name,Description)],by.x='gene_id',by.y='Name')
output=output%>%arrange(-rank)
output=output%>%dplyr::select(gene_name=Description,rank)



# output ranked list: 
write.table(output,file='../processed_data/160801/HCASMC_SF_vs_FBS.rnk',quote=F,row.names=F,sep='\t')


# save:
save(list = c('dds', 'res', 'resOrdered'), file = '../processed_data/160801/DESeq2_HCASMC_SF_vs_FBS.RData')
