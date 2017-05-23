#!/usr/bin/env Rscript
# bosh liu
# 2016/05/08
# durga
# run matrix eqtl with different set of PCs


# library
library("MatrixEQTL")
library('gap')

# paths:
figure_path='../figures/160516_matrix_eQTL/'


# command args: 
args=commandArgs(T)
SNP_file_name=args[1]
snps_location_file_name=args[2]
expression_file_name=args[3]
gene_location_file_name=args[4]
covariates_file_name=args[5]
output_dir=args[6]

# path:
SNP_file_name='../processed_data/160516_genotype/chr22.genotype.maf.txt'
snps_location_file_name='../processed_data/160516_genotype/chr22.genotype_loc.maf.txt'
expression_file_name='../processed_data/031_prepare_matrix_eQTL_expression/expression.txt'
gene_location_file_name='../processed_data/031_gen_gene_loc/gene_loc.txt'
covariates_file_name=""
output_dir='../processed_data/160516_genotype/'


# specify model:
useModel = modelLINEAR


# Output file name
output_file_name_cis = paste(output_dir,'cis.txt',sep="");
output_file_name_tra = paste(output_dir,'trans.txt',sep="");


# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 0;


# specify that the covariance matrix is a multiple of identity
errorCovariance = numeric()


# Distance for local gene-SNP pairs
cisDist = 1e6;


# read SNP file:
snps = SlicedData$new()
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name );



# read expression file: 
gene = SlicedData$new()
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
gene$LoadFile( expression_file_name );


# read covariates:
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
if(nchar(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}



## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
snps = snps, 
gene = gene, 
cvrt = cvrt,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel, 
errorCovariance = errorCovariance,
verbose = TRUE, 
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos, 
genepos = genepos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);



# results: 
head(me$cis$eqtls)
head(me$trans$eqtls)


# make Q-Q plot for cis eqtl: 
pdf(paste0(figure_path,'cis_eqtl_qqplot.pdf'))
qqunif(me$cis$eqtls$pvalue)
dev.off()


# now remove the three outliers: 
snps_mat=as.matrix(snps)
gene_mat=as.matrix(gene)
snps_mat2=snps_mat[,!colnames(snps_mat)%in%c('2135','2305','9070202')]
gene_mat2=gene_mat[,!colnames(gene_mat)%in%c('2135','2305','9070202')]


# create new sliced data objects: 
snps2=SlicedData$new()
snps2$CreateFromMatrix(snps_mat2)
gene2=SlicedData$new()
gene2$CreateFromMatrix(gene_mat2)


# run matrix eQTL:
me2 = Matrix_eQTL_main(
snps = snps2, 
gene = gene2, 
cvrt = cvrt,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel, 
errorCovariance = errorCovariance,
verbose = TRUE, 
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos, 
genepos = genepos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);


# the number of significant discoveries:
sum(me2$cis$eqtls$pvalue<0.05)
# 154424
sum(me$cis$eqtls$pvalue<0.05)
# 164573
sum(me2$cis$eqtls$FDR<0.05)  
# 1109
sum(me$cis$eqtls$FDR<0.05)
# 1746


# overlap qqplot for both 52 and 49 samples: 
png(paste0(figure_path,'52_vs_49_qqplot.png'))
qqplot(-log10(me2$cis$eqtls$pvalue),-log10(me$cis$eqtls$pvalue),main='p-value with 52 vs 49 samples',xlab='-log10(p-value): 49 samples',ylab='-log10(p-value): 52 samples')
abline(0,1)
dev.off()

# make histogram for p-value for 52 and 49 samples:
pdf(paste0(figure_path,'52_hist.pdf'))
hist(me$cis$eqtls$pvalue,breaks=100,main='p-value: 52 samples')
dev.off()

pdf(paste0(figure_path,'49_hist.pdf'))
hist(me2$cis$eqtls$pvalue,breaks=100,main='p-value: 49 samples')
dev.off()



# I am worried about the 1's in 49_hist.png:
# what are the 1's? 
idx=which(me2$cis$eqtls$pvalue>0.999)
snp=me$cis$eqtls$snps[idx[1]]
gen=me$cis$eqtls$gene[idx[1]]

pdf(paste0(figure_path,snp,'_',gen,'.pdf'))
plot(gene_mat2[(rownames(gene_mat2)==gen),],ylab='gene expression', main='49')
plot(gene_mat[(rownames(gene_mat)==gen),],ylab='gene_expression', main='52')
dev.off()
