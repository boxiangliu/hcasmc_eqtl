library(data.table)
hcasmc_gct='../processed_data/160715/combined.rpkm'
gtex_gct='/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/v6p/v6p_All_Tissues_gene_rpkm_FOR_QC_ONLY.gct'

hcasmc=fread(hcasmc_gct,header=T)
gtex=fread(gtex_gct,header=T)
merged=merge(gtex,hcasmc,by='Name')
merged[,Description.y:=NULL]
setnames(merged,'Description.x','Description')
fwrite(merged,file='../processed_data/160715/combined.gtex.hcasmc.rpkm',sep='\t',quote=F)