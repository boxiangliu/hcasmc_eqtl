#!/usr/bin/env Rscript
# boxiang liu
# durga
# find the top eqtls:


# command args:
args=commandArgs(T,T)
input_file=args$input
output_file=args$output
fdr=as.numeric(args$fdr)
# input_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160530/find_optimum_num_PEER_factors/fastqtl.pc3.peer8.padj.txt'
# output_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/compare_hcasmc_and_gtex/HCASMC.txt'
# fdr=0.05


# read input: 
hcasmc=fread(input_file,header=T)


# calculate the number of significant eQTLs: 
num_sig=sum(hcasmc$qvalues<=fdr)
message(num_sig,' eQTLs with qvalues <= ',fdr)


# find top eQTLs: 
top_eqtl=hcasmc[hcasmc$qvalues<=fdr,.(best_variant,gene_id)]

chr_pos=str_split_fixed(top_eqtl$best_variant,'_',5)[,c(1,2)]
top_eqtl$chr=chr_pos[,1]
top_eqtl$pos=chr_pos[,2]
top_eqtl$chr=str_replace(top_eqtl$chr,'chr','')
setcolorder(top_eqtl,c('chr','pos','gene_id','best_variant'))
top_eqtl$tissue='HCASMC'


# write output:
write.table(top_eqtl,file=output_file,quote=F,col.names=F,row.names=F,sep='\t')


