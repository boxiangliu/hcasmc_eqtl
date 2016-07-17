#!/usr/bin/env Rscript
# boxiang liu
# durga
# get top GWAS hits


# command args:
args=commandArgs(T,T)
input_file=args$input
output_file=args$output
alpha=as.numeric(args$alpha)
# input_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/gwas/cad.add.160614.website.txt'
# output_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/compare_hcasmc_and_gtex/GWAS.txt'
# alpha=1e-7

# read input:
input=fread(input_file)



# calculate the number of significant variants: 
num_sig=sum(input$p_dgc<=alpha)
message(num_sig,' GWAS hits with p-value <= ', alpha)


# subset to significant hits: 
output=input[input$p_dgc<=alpha,]
output$variant_id=with(output,paste(chr,bp_hg19,effect_allele,noneffect_allele,'b37',sep='_'))
output=output[,.(chr,bp_hg19,markername,variant_id)]
output$annotation='GWAS'


# write output:
write.table(output,file=output_file,quote=F,sep='\t',col.names=F,row.names=F)