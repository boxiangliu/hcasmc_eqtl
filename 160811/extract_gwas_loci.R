# extract GWAS hits from the supplmentary excel in the paper Nikpay 2015 Nature Genetics
library(XLConnect)

# read gwas loci:
gwas_loci=readWorksheet(loadWorkbook("/srv/persistent/bliu2/HCASMC_eQTL/data/gwas/nikpay_2015_ng.xlsx"),sheet=6,startRow=3,endRow=59,startCol=1,endCol=9,check.names=F,header=T,)


# set column order:
setcolorder(gwas_loci,c(2,4,1,3,seq(5,9)))


# write gwas loci:
write.table(gwas_loci,'../processed_data/160811/gwas_loci.cad.all.FWER.txt',sep='\t',row.names=F,quote=F)