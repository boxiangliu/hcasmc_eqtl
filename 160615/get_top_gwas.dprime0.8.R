#!/usr/bin/env Rscript
# boxiang liu
# durga
# get top GWAS hits


# command args:
args=commandArgs(T,T)
input_file=args$input
output_file=args$output
lead_snp_file=args$lead
# input_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/gwas/cad.add.160614.website.txt'
# lead_snp_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/gwas/CARDIOGRAMplusC4DleadSNPsplusSNPsLD.Dprime0.8'
# output_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/compare_hcasmc_and_gtex/GWAS.txt'


# read input:
input=fread(input_file)
lead_snp=fread(lead_snp_file,header=T)


# subset to lead SNPs:
output=input[markername%in%lead_snp$Proxy,]
output$variant_id=with(output,paste(chr,bp_hg19,effect_allele,noneffect_allele,'b37',sep='_'))
output=output[,.(chr,bp_hg19,markername,variant_id)]
output$annotation='GWAS'


# write output:
write.table(output,file=output_file,quote=F,sep='\t',col.names=F,row.names=F)