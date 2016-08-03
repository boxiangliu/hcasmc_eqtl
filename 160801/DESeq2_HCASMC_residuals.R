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


# remove the singleton 9090701: 
count_wide=count_wide%>%dplyr::select(-`9090701_sf`)


# coerce to data.frame: 
count_matrix = count_wide %>% dplyr::select(-Name,-Description) %>% data.frame(check.names=F)
rownames(count_matrix) = count_wide %>% select(Name) %>% unlist()


# filter out genes with less than 1 read count: 
count_matrix2=count_matrix[rowSums(count_matrix)>=1,]

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
  batch = batch
)


########################
# create DESeq2 dataset
########################
dds = DESeqDataSetFromMatrix(countData = count_matrix2,
                             colData = colData,
                             design = ~ batch)


##############
# run DESeq2 
##############
dds = DESeq(dds, parallel = TRUE)
# sum(resOrdered$padj<1e-5,na.rm=T)

fitted.common.scale = t(t(assays(dds)[["mu"]])/sizeFactors(dds))
residuals=counts(dds, normalized=TRUE)-fitted.common.scale


# format the output: 
col_data=merge(x=data.frame(Name=rownames(residuals)),
	  y=count_wide[,.(Name,Description)],
	  by='Name')
output=data.frame(Name=rownames(residuals),residuals,check.names=F)
output=merge(col_data,output,by='Name')


# switch the name and the description column: 
setnames(output,c('Name','Description'),c('Description','Name'))
setcolorder(output,colnames(output)[c(2,1,3:ncol(output))])


# output residuals: 
write.table(output,file='../processed_data/160801/HCASMC_SF_vs_FBS.txt',quote=F,row.names=F,sep='\t')


# output phenotype labels:
condition
nsamples=ncol(residuals)
nclasses=2
if (file.exists('../processed_data/160801/HCASMC_SF_vs_FBS.cls')) file.remove('../processed_data/160801/HCASMC_SF_vs_FBS.cls')
write(paste(nsamples,nclasses,'1',sep='\t'),file='../processed_data/160801/HCASMC_SF_vs_FBS.cls',sep='\t',append=T)
write("#\tfbs\tsf",file='../processed_data/160801/HCASMC_SF_vs_FBS.cls',sep='\t',append=T)
write(condition,ncolumns=nsamples,file='../processed_data/160801/HCASMC_SF_vs_FBS.cls',sep='\t',append=T)
