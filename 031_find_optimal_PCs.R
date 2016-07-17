#!/usr/bin/env Rscript
# bosh liu
# 2016/05/08
# durga
# run matrix eqtl with different set of PCs


# library
library("MatrixEQTL")
library('gap')


# path: 
input_dir='../processed_data/031_find_optimal_PCs/'
output_dir='../processed_data/031_find_optimal_PCs/'

# specify model:
useModel = modelLINEAR


# specify input file names: 
SNP_file_name = paste(input_dir, "chr20.genotype.txt", sep="");
expression_file_name = paste(input_dir, "expression.txt", sep="");
covariates_file_name = paste(input_dir, "covariates.txt", sep=""); 
snps_location_file_name = paste(input_dir, "chr20.genotype_loc.txt", sep="");
gene_location_file_name = paste(input_dir, "chr20.gene_loc.txt", sep="");

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
if(length(covariates_file_name)>0) {
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
pdf('1.pdf')
qqunif(me$cis$eqtls$pvalue)
dev.off()


# there are a few off diagonal points, and 
# several things needs to be checked: 
# 1) normality of residuals (for significant points) 
# 2) plot the genotype versus the expression for the highly significant points 
gene_matrix=as.matrix(gene)
snp_matrix=as.matrix(snps)
pdf('4.pdf')
for (i in 1:10){
	# get index for most significatn SNP-gene pair:
	idx=i

	# get expression for the most significant association:
	gene_name=me$cis$eqtls[idx,'gene']
	gene_row=gene_matrix[gene_name,]


	# get genotype for the most significant association: 
	snp_name=me$cis$eqtls[idx,'snps']
	snp_row=snp_matrix[snp_name,]


	# plot expression vs genotype: 
	plot(snp_row, gene_row,main=i)
	# trend does not look significant.. 
}
dev.off()


idx=10000

# get expression for the most significant association:
gene_name=me$cis$eqtls[idx,'gene']
gene_row=gene_matrix[gene_name,]

head(me$cis$eqtls)



# get genotype for the most significant association: 
snp_name=me$cis$eqtls[idx,'snps']
snp_row=snp_matrix[snp_name,]


# plot expression vs genotype: 
plot(snp_row, gene_row,main=idx)

# get index for second most significatn SNP-gene pair:
idx=2

# get expression for the most significant association:
gene_name=me$cis$eqtls[idx,'gene']
gene_matrix=as.matrix(gene)
gene_row=gene_matrix[gene_name,]


# get genotype for the most significant association: 
snp_name=me$cis$eqtls[idx,'snps']
snp_matrix=as.matrix(snps)
snp_row=snp_matrix[snp_name,]


# plot expression vs genotype: 
plot(snp_row, gene_row)

# make histogram for cis eqtl: 
pdf('2.pdf')
hist(me$cis$eqtls$pvalue)
dev.off()


# make histogram for trans eqtl: 
pdf('3.pdf')
hist(me$trans$eqtls$pvalue)
dev.off()


