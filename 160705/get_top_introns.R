# subset top introns based on the number of supporting reads: 

# command args: 
args=commandArgs(T,T)
ratio_file=args$ratio
count_file=args$count
output_file=args$output
ratio_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160627/leafcutter.norm.tsv'
count_file='/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/leafcutter/leafcutter_perind_numers.counts.gz'
output_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160705/leafcutter.norm.top10000.tsv'


# read input: 
ratio=read.table(ratio_file,row.names=1,header=1,check.names=F)
count=read.table(count_file,check.names=F)


# remove sex chromosomes: 
count=count[!str_detect(row.names(count),'chrY|chrX'),]
ratio=ratio[!str_detect(row.names(ratio),'chrY|chrX'),]
stopifnot(dim(ratio)==dim(count))


# reorder rows:
idx=match(row.names(ratio),row.names(count))
count=count[idx,]
stopifnot(row.names(ratio)==row.names(count))


# subset top introns
row_sums=rowSums(count)
sum_rank=rank(-row_sums,ties.method='min')
ratio_top=ratio[sum_rank%in%seq(1,10000),]
ratio_top=data.table(Name=row.names(ratio_top),ratio_top)


# write output:
write.table(ratio_top,file=output_file,quote=F,sep='\t',row.names=F)
