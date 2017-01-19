# Libraries: 
library(data.table)
library(dplyr)
library(dtplyr)
library(cowplot)
library(stringr)

# Variables:
ecaviar_variants_file='../processed_data/mpra/forNathan/eCAVIAR_variants.bed'
naive_overlap_file='../processed_data/mpra/forNathan/naive_overlap_variants.bed'
ecaviar_variants_outfile='../processed_data/mpra/forNathan/eCAVIAR_variants.clpp5e-3.bed'
naive_overlap_outfile1='../processed_data/mpra/forNathan/naive_overlap_variants.p1e-6.rank50.bed'
naive_overlap_outfile2='../processed_data/mpra/forNathan/naive_overlap_variants.p1e-6.rank40.bed'

# Read eCAVIAR variants: 
ecaviar_variants=fread(ecaviar_variants_file)


# Filter to variants in both GWAS and eQTL causal set: 
ecaviar_variants%>%filter(gwasSet==1 & eqtlSet==1)%>%nrow()
# 2534 

# Filter to variants above CLPP cutoff: 
ecaviar_variants%>%filter(clpp>0.01)%>%nrow() # 130
ecaviar_variants%>%filter(clpp>0.005)%>%nrow() # 252
ecaviar_variants%>%filter(clpp>0.001)%>%nrow() # 821

# Output variants above CLPP 0.005: 
write.table(ecaviar_variants%>%filter(clpp>0.005),ecaviar_variants_outfile,sep='\t',row.names=F,col.names=T,quote=F)


# Read naive overlap variants: 
naive_overlap=fread(naive_overlap_file)


# Filter to variants above p-value cutoff: 
naive_overlap%>%filter(pval<1e-6,rank<=30,str_detect(gene_id,'ENSG00000213445'))
naive_overlap%>%filter(pval<1e-6,rank<=30,str_detect(gene_id,'ENSG00000118526')) 
naive_overlap%>%filter(pval<1e-6,rank<=30,str_detect(gene_id,'ENSG00000166949')) 
naive_overlap%>%filter(pval<1e-6,rank<=30)%>%nrow() # 115
naive_overlap%>%filter(pval<1e-6,rank<=40)%>%nrow() # 148
naive_overlap%>%filter(pval<1e-6,rank<=50)%>%nrow() # 253


# Output to variants below p-value 1e-6 and are top 50 variant: 
write.table(naive_overlap%>%filter(pval<1e-6,rank<=50),naive_overlap_outfile1,sep='\t',row.names=F,col.names=T,quote=F)
write.table(naive_overlap%>%filter(pval<1e-6,rank<=40),naive_overlap_outfile2,sep='\t',row.names=F,col.names=T,quote=F)