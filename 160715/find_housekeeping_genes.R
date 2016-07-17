#!/usr/bin/env Rscript
# boxiang liu
# durga
# differential expression 

# library:
source('/srv/persistent/bliu2/HCASMC_eQTL/scripts/utils.R')
library("BiocParallel")
register(MulticoreParam(40))
library('DESeq2')


# read input:
hcasmc_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/combined.count'
hcasmc=read.table(hcasmc_file,header=T,check.names=F)
res=decompose_gct(hcasmc)



# inititilize: 
master_count=res$count
master_col_data=data.frame(res$col_data,tissue='HCASMC')
master_row_data=res$row_data


gtex_files=list.files('/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/v6p/subsampling',pattern='*.52.count',recursive=T,full.name=T)
for (gtex_file in gtex_files){
	tissue=str_replace(basename(gtex_file),'_subsample.52.count','')
	message(tissue)
	count=read.table(gtex_file,header=T,check.names=F)

	res=decompose_gct(count)
	count=res$count
	col_data=res$col_data
	row_data=res$row_data
	col_data$tissue=tissue


	# sanity check:
	stopifnot(setequal(row_data$Name,master_row_data$Name))
	stopifnot(length(unique(row_data$Name))==length(master_row_data$Name))


	# reorder the rows of count: 
	idx=match(rownames(master_count),rownames(count))
	count=count[idx,]
	row_data=row_data[idx,]
	stopifnot(row_data$Name==master_row_data$Name)
	stopifnot(row.names(count)==row.names(master_count))


	# append to master: 
	master_col_data=rbind(master_col_data,col_data)
	master_count=data.frame(master_count,count,check.names=F)
}


# subset the master_count: 

master_count

# create DESeq dataset: 
dds=DESeqDataSetFromMatrix(countData = as.matrix(master_count),colData=master_col_data,design=~tissue)
dds=dds[rowSums(counts(dds))>=50,]


# run DESeq:
dds=DESeq(dds,test='LRT',reduced=~1,parallel=T)
res=results(dds)
resOrdered=res[order(res$padj),]


# save: 
save(dds,res,resOrdered,file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/find_housekeeping_genes.RData')

