# command args: 
args=commandArgs(T,T)
input=args$input
output=args$output
input='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160627/leafcutter.norm.bed.gz'
output='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160705/leafcutter.norm.perm.bed'



# read input: 
leafcutter=fread(sprintf('zcat %s',input),header=T)


# save column name order: 
col_names=colnames(leafcutter)


# permute columns: 
set.seed(1)
idx=sample(5:ncol(leafcutter),ncol(leafcutter)-4)
leafcutter_perm=leafcutter[,c(1:4,idx),with=F]


# change column names back:
setnames(leafcutter_perm,col_names)

# write output:
write.table(leafcutter_perm,file=output,quote=F,row.names=F,sep='\t')


