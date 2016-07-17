library(stringr)
library(cowplot)
library(data.table)
#------ 160715 -----
#' decompose a gct count file:
#' a gct count files' first and second columns are ENSEMBL ID and gene name
#' @return res (list) with three elements: count, col_data and row_data 
decompose_gct=function(gct){
	gct=as.data.frame(gct)
	row_data=gct[,1:2]
	col_data=data.frame(sample=colnames(gct[,-c(1,2)]))
	count=gct[,-c(1,2)]
	rownames(count)=gct[,1]
	res=list(count=count,col_data=col_data,row_data=row_data)
	return(res)
}


#------ 160527 -----
#' inverse rank normalization
getInverseNormal <- function(x){
	x <- as.numeric(x)
	xx2<-qnorm((rank(x,na.last = "keep") - 0.5) / sum(!is.na(x)))
	return(xx2)
}

#------ 30_compare_9052004_samples.R ------
#' make pairs upper panel show correlation:
#' taken from help file of pairs()
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    # par(usr = c(0, 1, 0, 1))
    par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)

    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}


#------ convert_imput2_output_to_genotype.R -----
#' convert impute2 output to genotype:
#' @param impute2 (data.table(i,j)) impute2 output
#' @return genotypes (matrix(i,(j-5)/3)) genotypes
convImpute2ToGeno=function(impute2){
	n_samples=(ncol(impute2)-5)/3
	stopifnot(as.integer(n_samples)==n_samples)
	n_snps=nrow(impute2)
	genotypes=matrix(ncol=n_samples,nrow=n_snps)
	colnames(genotypes)=1:n_samples
	for (j in seq(1,n_samples)){
		temp=impute2[,((j-1)*3+1+5):(j*3+5),with=F]
		genotypes[,j]=max.col(temp)-1
	}
	return(genotypes)
}


#' report error if multi-allelic sites are found. 
#' @param genotypes (data.table) genotypes with samples as columns and variants as rows
#' @return reports error and aborts if multi-allelic loci are found. 
detectMultiAllelicLoci = function(genotypes){
	genotypes_long = melt(genotypes, id.vars = c('CHROM','POS'), variable.name = 'ID', value.name = 'genotype')
	genotype_list = names(table(genotypes_long$genotype))
	if (!all(genotype_list %in% c('0/0','0/1','1/1','0|0','0|1','1|0','1|1'))) {
		stop('multi-allelic loci found!')
	} else {
		message('no multi-allelic loci.')
	}
}


#' remove multi allelic loci (e.g. 0/2, 1/2)
#' @param genotypes (data.table) genotypes with samples as columns and variants as rows
#' @return (data.table) genotypes with only the following genotypes: 0/0, 0/1, 1/1, 0|0, 0|1, 1|0, 1|1. 
removeMultiAllelicLoci = function(genotypes){
	genotypes_long = melt(genotypes, id.vars = c('CHROM','POS'), variable.name = 'sample', value.name = 'genotype')
	genotype_counts = table(genotypes_long$genotype)
	all_genotypes = names(genotype_counts)
	multiallelic_genotypes = all_genotypes[!all_genotypes %in% c('0/0','0/1','1/1','0|0','0|1','1|0','1|1')]

	for (s in multiallelic_genotypes){
		multi = unique(which(genotypes == s, arr.ind = T)[,'row'])
		if (length(multi) == 0){
			message(sprintf("%s not found in the remaining loci",s))
		} else {
			genotypes = genotypes[-multi,]
		}
	}
	return(genotypes)
}




#' convert genotypes to dosages
#' @param genotypes (data.table) genotypes with samples as columns and variants as rows; should only have biallalic variants (e.g. 0/0, 0/1, 1/1); the first and second columns should be 'CHROM' and 'POS'. 
#' @return (numeric matrix) dosages correspond to the genotypes; 'CHROM' and 'POS' column are turned into rownames. 
genotypeToDosage = function(genotypes){
	stopifnot('data.table' %in% class(genotypes))
	for (s in c('0/2','1/2','2/2','0|2','2|0','1|2','2|1','2|2')) {stopifnot(length(which(genotypes == s)) == 0)}
	# genotype_long = melt(data = genotypes, id.vars = c('CHROM','POS'))
	dosages = matrix(0, nrow = nrow(genotypes), ncol = ncol(genotypes))
	colnames(dosages) = colnames(genotypes)
	rownames(dosages) = genotypes[,paste(CHROM, POS, sep = '_')]
	for (s in c('0/0','0/1','1/1','0|0','0|1','1|0','1|1')){
		dosage = as.integer(strsplit(s, split = "/|\\|")[[1]])
		dosage = sum(dosage)
		idx = which(genotypes == s)
		dosages[idx] = dosage
	}
	dosages = dosages[,!(colnames(dosages) %in% c('CHROM','POS'))]
	return(dosages)
}



#' remove missing genotypes
#' @param genotypes (data.table) genotypes with samples as columns and variants as rows
#' @return genotypes but rows with missing values (i.e. ./.) are removed 
removeMissingGenotypes = function(genotypes){
	missing = unique(which(genotypes == "./.", arr.ind=T)[,'row'])
	stopifnot(length(unique(missing)) == length(missing))
	if (length(missing) == 0) {
		message('no missing genotypes')
	} else {
		genotypes = genotypes[-missing,]
	} # remove rows with missing genotypes
	return(genotypes)
}

#' read sample list
#' @param filename (character) a file with one line per sample
#' @return (character vector) of samples
readSampleList = function(filename){
	sample_list = read.table(filename, header = FALSE, stringsAsFactors = FALSE)
	sample_list = unlist(sample_list)
	stopifnot(class(sample_list) == 'character')
	return(sample_list)
}

#' calculate reference ratio 
#' @param sample (data.table) table with columns refCount, refCount, totalCount, other columns will be ignored. 
#' @return a data.table with an addition column refRatio
calcRefRatio = function(sample){
	stopifnot('data.table' %in% class(sample))
	sample[, refRatio := refCount/totalCount]
	return(sample[])
}


#' @param assays (data.table) table with columns ID (variant ID), type (alt, ref or totalCount), and count. 
#' @param sample_name the name of the sample
#' @return a wide data.table with columns ID, refCount, altCount, totalCount 
getSample = function(assays, sample_name){
	stopifnot('data.table' %in% class(assays))
	sample = assays[sample == sample_name, ]
	sample = data.table::dcast(sample, ID ~ type, value.var = 'count')
	return(sample)
}

#' @param assays table with columns ID (variant ID), type (alt, ref or totalCount), and count. 
#' @param countType one of totalCount, altCount or refCount
#' @return a wide data.table with an ID column, and column for each sample
getCount = function(assays, countType = c('totalCount', 'altCount','refCount')){
	countType = match.arg(countType)
	stopifnot('data.table' %in% class(assays))
	count = assays[type == countType, ]
	count = data.table::dcast(count, ID ~ sample, value.var = 'count')
	return(count)
}

