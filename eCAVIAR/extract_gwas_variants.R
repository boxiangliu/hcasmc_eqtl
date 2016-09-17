# extract GWAS hits from the supplmentary excel in the paper Nikpay 2015 Nature Genetics
library(XLConnect)

# command args: 
args=commandArgs(T)
in_file=args[1]
out_file=args[2]
# in_file="/srv/persistent/bliu2/HCASMC_eQTL/data/gwas/nikpay_2015_ng.xlsx"
# out_file='../processed_data/eCAVIAR/gwas_loci.cad.all.genomewide_fdr_merged.txt'


# read gwas loci:
gwas_fdr_hits=readWorksheet(loadWorkbook(in_file),sheet=7,startRow=2,endRow=204,startCol=1,endCol=3,check.names=F,header=T)%>%select(markername,chr,pos=bp_hg19)
gwas_genowide_hits=readWorksheet(loadWorkbook(in_file),sheet=6,startRow=3,endRow=59,startCol=2,endCol=4,check.names=F,header=T)%>%select(markername=SNP,chr=CHR,pos=`Base-pair Position`)


# merge two lists: 
merged=rbind(gwas_fdr_hits,gwas_genowide_hits)

# take unique hits: 
uniq=unique(merged)

# output the merged unique hits:
write.table(uniq,out_file,sep='\t',row.names=F,quote=F)