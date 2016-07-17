#!/usr/bin/env Rscript
# boxiang liu
# durga
# find the top eqtls


# command args:
args=commandArgs(T,T)
gtex_file=args$input
output_dir=args$outdir
num_sig=as.numeric(args$num)

# gtex_file='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6_fastQTL_FOR_QC_ONLY/Uterus_Analysis.v6.FOR_QC_ONLY.egenes.txt.gz'
# output_dir='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/compare_hcasmc_and_gtex/'
# num_sig=290



# read input: 
message(gtex_file)
options(warn=-1)
gtex=fread(sprintf('zcat %s',gtex_file),header=T)
options(warn=0)

# get tissue name: 
tissue=basename(gtex_file)%>%str_replace('_Analysis.v6.FOR_QC_ONLY.egenes.txt.gz','')


# remove sex chromosome: 
gtex=gtex[gtex$gene_chr!='X',]


# get top eQTLs:
top_eqtl=gtex[gtex$qval%in%sort(gtex$qval)[1:num_sig],.(variant_id,gene_id)]
chr_pos=str_split_fixed(top_eqtl$variant_id,'_',5)[,c(1,2)]
top_eqtl$chr=chr_pos[,1]
top_eqtl$pos=chr_pos[,2]
setcolorder(top_eqtl,c('chr','pos','gene_id','variant_id'))
top_eqtl$tissue=tissue


# write output:
output_file=paste0(output_dir,'/',tissue,'.txt')
write.table(top_eqtl,file=output_file,quote=F,col.names=F,row.names=F,sep='\t')

