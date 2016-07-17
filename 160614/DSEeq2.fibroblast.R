#!/usr/bin/env Rscript
# boxiang liu
# durga
# differential expression 

# library:
library("BiocParallel")
register(MulticoreParam(40))
library('DESeq2')


# command args: 
args=commandArgs(T,T)
gtex_file=args$gtex
gtex_col_data_file=args$gtex_col_data
hcasmc=args$hcasmc
# gtex_file='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6p_All_Tissues_read_counts_FOR_QC_ONLY.gct'
# gtex_col_data_file='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/sample_annotations/GTEx_Analysis_2015-01-12_Annotations_SampleAttributesDS.txt'
# hcasmc_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160614/rnaseqc.hcasmc_eqtl.reads.gct'


# read input
gtex=fread(gtex_file,skip=2,header=T)
hcasmc=fread(hcasmc_file,skip=2,header=T)
gtex_col_data=fread(gtex_col_data_file,header=T)


# sanity check:
stopifnot(setequal(gtex$Name,hcasmc$Name))
stopifnot(length(unique(gtex$Name))==length(gtex$Name))


# subset to one tissue:
tissue_name='Cells - Transformed fibroblasts'
tissue_idx=gtex_col_data[SMTSD==tissue_name,SAMPID]
tissue=gtex[,c(1,2,which(colnames(gtex)%in%tissue_idx)),with=F]


# create count matrix:
merged=merge(tissue,hcasmc[,-2,with=F],by='Name')
idx=match(tissue$Name,merged$Name)
merged=merged[idx,]
count=as.matrix(merged[,-c(1,2),with=F])
rownames(count)=merged$Name


# create column data:
coldata=data.frame(tissue=c(rep(tissue_name,ncol(tissue)-2),rep('HCASMC',ncol(hcasmc)-2)))


# create DESeq dataset: 
dds=DESeqDataSetFromMatrix(countData = count,colData=coldata,design=~tissue)
dds=dds[rowSums(counts(dds))>1,]


# run DESeq:
dds=DESeq(dds,parallel=T)
res=results(dds)
resOrdered=res[order(res$padj),]


# number of significant genes:
FDR=0.05
num_sig=sum(res$padj<FDR, na.rm=TRUE)
message(tissue_name,' vs HCASMC: ',num_sig,' DE genes at FDR=',FDR)



# make MA plot:
pdf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160614/ma.pdf')
plotMA(res, main="DESeq2", ylim=c(-2,2))
dev.off()


# plot fibroblast specific genes: 
# ENSG00000196154.7 (fibroblasts specific gene):
gene_id='ENSG00000196154.7'
pdf(sprintf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160614/%s.pdf',gene_id))
plotCounts(dds, gene=gene_id, intgroup="tissue", main='S100A4 (fibroblast specific gene)')
dev.off()


# Collagen 
gene_id='ENSG00000164692.13'
pdf(sprintf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160614/%s.pdf',gene_id))
plotCounts(dds, gene=gene_id, intgroup="tissue", main='COL1A2')
dev.off()

gene_id='ENSG00000168542.8'
pdf(sprintf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160614/%s.pdf',gene_id))
plotCounts(dds, gene=gene_id, intgroup="tissue", main='COL3A1')
dev.off()


# thymic stromal lymphopoietin (TSLP, ENSG00000145777.10):
# induces the release of T cell-attracting chemokines
gene_id='ENSG00000145777.10'
pdf(sprintf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160614/%s.pdf',gene_id))
plotCounts(dds, gene=gene_id, intgroup="tissue", main='TSLP')
dev.off()


# plot vimentin (ENSG00000026025.9)
gene_id='ENSG00000026025.9'
pdf(sprintf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160614/%s.pdf',gene_id))
plotCounts(dds, gene=gene_id, intgroup="tissue", main='Vimentin')
dev.off()


# plot smooth muscle cell specific genes: 
# ENSG00000122786.15, caldesmon 1
gene_id='ENSG00000122786.15'
gene_name=gtex[Name==gene_id,Description]
pdf(sprintf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160614/%s.pdf',gene_id))
plotCounts(dds, gene=gene_id, intgroup="tissue", main=gene_name)
dev.off()


# ENSG00000133026.8 (myosin,MYH10):
gene_id='ENSG00000133026.8'
gene_name=gtex[Name==gene_id,Description]
pdf(sprintf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160614/%s.pdf',gene_id))
plotCounts(dds, gene=gene_id, intgroup="tissue", main=gene_name)
dev.off()


# ENSG00000167460.10 (tropomyosin,TPM4):
gene_id='ENSG00000167460.10'
gene_name=gtex[Name==gene_id,Description]
pdf(sprintf('/srv/persistent/bliu2/HCASMC_eQTL/figures/160614/%s.pdf',gene_id))
plotCounts(dds, gene=gene_id, intgroup="tissue", main=gene_name)
dev.off()