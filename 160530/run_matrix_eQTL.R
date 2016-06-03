#!/usr/bin/env Rscript
# bosh liu
# 2016/05/08
# durga
# run matrix eqtl with different set of PCs


# library
library("MatrixEQTL")
library('gap')


# command args: 
args=commandArgs(T)
SNP_file_name=args[1]
snps_location_file_name=args[2]
expression_file_name=args[3]
gene_location_file_name=args[4]
covariates_file_name=args[5]
output_dir=args[6]

# path:
# SNP_file_name='../processed_data/160516_genotype/chr22.genotype.maf.txt'
# snps_location_file_name='../processed_data/160516_genotype/chr22.genotype_loc.maf.txt'
# expression_file_name='../processed_data/031_prepare_matrix_eQTL_expression/expression.txt'
# gene_location_file_name='../processed_data/031_gen_gene_loc/gene_loc.txt'
# covariates_file_name=""
# output_dir='../processed_data/160516_genotype/'


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
output_file_name = output_file_name_tra,
pvOutputThreshold = pvOutputThreshold_tra,
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



# make Q-Q plot for cis eqtl: 
pdf(paste0(output_dir,'/cis_eqtl_qqplot.pdf'))
qqunif(me$cis$eqtls$pvalue)
dev.off()


