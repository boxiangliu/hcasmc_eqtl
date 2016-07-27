# find housekeeping genes with mean and variance analysis: 

# libraries:
source('/srv/persistent/bliu2/HCASMC_eQTL/scripts/utils.R')


# read input:
hcasmc_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160519_rpkm/combined.rpkm'
hcasmc=read.table(hcasmc_file,header=T,check.names=F)
hcasmc_sub=subsample_gct(hcasmc,10)
res=decompose_gct(hcasmc_sub)


# inititilize: 
master_count=res$count
master_col_data=data.frame(res$col_data,tissue='HCASMC')
master_row_data=res$row_data


gtex_files=list.files('/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/v6p/subsampling',pattern='*.10.rpkm',recursive=T,full.name=T)
for (gtex_file in gtex_files){
	tissue=str_replace(basename(gtex_file),'_subsample.10.rpkm','')
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


# merge count table and col_data to facilitate melting:
master_count_t=data.frame(t(master_count))
master_count_t$sample=rownames(master_count_t)
merged=merge(master_col_data,master_count_t,by='sample')


# melt merged data to facilitate filtering:
melted=melt(merged,id=c('sample','tissue'),variable.name='gene_id',value.name='rpkm')
melted=melted%>%group_by(gene_id,tissue)%>%mutate(median=median(rpkm))
melted=melted%>%dplyr::select(tissue,gene_id,median)%>%unique()
melted=as.data.table(melted)


# calculate log median rpkm and sd: 
melted[,logmedian:=log2(median+1)]
melted[,sd:=sd(logmedian),by='gene_id']


# apply sd filter, expression filter, and log fold change filter:
filtered=melted[sd<1,]
filtered[,min:=min(median),by='gene_id']
filtered=filtered[min>0,]
filtered[,avg:=mean(logmedian),by='gene_id']
filtered[,maxdiff:=max(logmedian)-avg,by='gene_id']
filtered=filtered[maxdiff<1,]


# get and save housekeeping gene names:
hk_genes=unique(filtered$gene_id)
hk_genes=hcasmc[hcasmc$Name%in%hk_genes,c('Name','Description')]
write.table(hk_genes,'/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/hk_genes.txt',quote=F,sep='\t',row.names=F)


# plot a housekeeping gene: 
hk_gene_exp=master_count[rownames(master_count)==as.character(hk_genes[1,1]),]
stopifnot(master_col_data$sample==colnames(hk_gene_exp))
pdf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/read_counts_for_hk_gene.pdf')
to_plot=data.frame(tissue=master_col_data$tissue,rpkm=unlist(hk_gene_exp))
ggplot(to_plot,aes(tissue,rpkm))+geom_boxplot()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()


# plot a random non housekeeping gene: 
non_hk_gene_exp=master_count[which(!rownames(master_count)%in%hk_genes[,1])[2],]
pdf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/read_counts_for_non_hk_gene.pdf')
to_plot=data.frame(tissue=master_col_data$tissue,rpkm=unlist(non_hk_gene_exp))
ggplot(to_plot,aes(tissue,rpkm))+geom_boxplot()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()