# library: 
library(dplyr)

# command args: 
args=commandArgs(T)
in_file=args[1]
out_dir=args[2]
# in_file='../processed_data/eCAVIAR/gwas_loci.cad.all.genomewide_fdr_merged.txt'
# out_dir='../processed_data/eCAVIAR/cad_gwas_loci/'


# read gwas data: 
gwas=fread('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/gwas/cad.add.160614.website.txt')
gwas=gwas%>%arrange(chr,bp_hg19)%>%as.data.table()
gwas=gwas[,index:=seq(nrow(gwas))]


# read top variants: 
top=fread(in_file)
for (i in seq(nrow(top))){
	variant=top[i,]
	message(variant$markername)
	idx=gwas[(chr==variant$chr)&(bp_hg19==variant$pos),index]
	out_file=paste0(variant$chr,"_",variant$pos,"_",variant$markername,'.txt')
	write.table(gwas[(idx-50):(idx+50),],file=paste0(out_dir,out_file),sep='\t',quote=F,row.names=F,col.names=T)
}