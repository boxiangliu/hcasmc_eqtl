# bosh liu
# write column data (sample, tissue) to a text file: 

# read col_data:
load('../processed_data/160715/find_housekeeping_genes.RData')
col_data=as.data.frame(colData(dds)[,c('sample','tissue')])


# write col_data:
write.table(col_data,file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160715/col_data.txt',quote=F,row.names=F,sep="\t")
