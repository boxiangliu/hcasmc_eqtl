library(XLConnect)
library(data.table)
library(DESeq2)
source('/srv/persistent/bliu2/rpe/scripts/utils/genome_annotation.R')
source('rebuttal/utils.R')

count_fn = '../data/rnaseq2/read_count/rnaseqc/rnaseqc.hcasmc_eqtl.reads.gct'

read_count = function(count_fn){
	count = fread(count_fn,header=TRUE)
	setDF(count)
	rownames(count) = count$Name
	count$Name = NULL
	count$Description = NULL
	return(count)
}

covariate = read_known_covariate(covariate_fn)

count = read_count(count_fn)
covariate = covariate[match(colnames(count),covariate$RNA),]

dds = DESeqDataSetFromMatrix(countData = count, colData = covariate, design = ~ Genomic_Ethnicity+Sex+Age)
dds = DESeq(dds)
age_res = as.data.frame(results(dds))
age_res[age_res$padj<0.05&!is.na(age_res$padj),]


dds = DESeqDataSetFromMatrix(countData = count, colData = covariate, design = ~ Genomic_Ethnicity+Age+Sex)
dds = DESeq(dds)
sex_res = as.data.frame(results(dds))
sig_sex_res = sex_res[sex_res$padj<0.05&!is.na(sex_res$padj),]
nrow(sig_sex_res) # 90
setorder(sig_sex_res,padj)

gene_annotation = read_gencode()
sig_sex_res$gene_id = rownames(sig_sex_res)
sig_sex_res=merge(sig_sex_res,gene_annotation[,list(gene_id,gene_name,chr)],,by='gene_id',sort=FALSE)
sum(sig_sex_res$chr == 'chrX') # 15
sum(sig_sex_res$chr == 'chrY') # 21