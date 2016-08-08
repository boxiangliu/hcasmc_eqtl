# read in HCASMC:
hcasmc=fread('zcat /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/hcasmc.eqtl.pc4.peer8.txt.gz')
setnames(hcasmc,c('pheno','geno','dist','pval','beta','se'))


# reformat the geno column: 
geno=str_replace(hcasmc$geno,'chr','')
geno=paste(geno,'b37',sep='_')
hcasmc$geno=geno

# write output: 
write.table(hcasmc,'/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/hcasmc.eqtl.pc4.peer8.b37.txt',sep='\t',row.names=F,col.names=F,quote=F)