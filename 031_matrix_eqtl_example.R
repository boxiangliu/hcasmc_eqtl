#!/usr/bin/env Rscript
# bosh liu
# 2016/05/08
# durga
# some matrix eQTL sample code


# library
library("MatrixEQTL")


# load example data: 
base.dir = find.package("MatrixEQTL")


# specify model:
useModel = modelLINEAR


# specify input file names: 
SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");
expression_file_name = paste(base.dir, "/data/GE.txt", sep="");
covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep=""); 
snps_location_file_name = paste(base.dir, "/data/snpsloc.txt", sep="");
gene_location_file_name = paste(base.dir, "/data/geneloc.txt", sep="");

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();


# Only associations significant at this level will be saved
pvOutputThreshold_cis = 2e-2;
pvOutputThreshold_tra = 1e-2;


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

unlink(output_file_name_tra);
unlink(output_file_name_cis);


# results: 
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)
