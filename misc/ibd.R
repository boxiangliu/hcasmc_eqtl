#!/usr/bin/env Rscript
# bosh liu
# 2016/04/19
# durga
# make IBD2* ratio plot 

source('utils.R')

genotype_filename = "../data/joint/recalibrated_variants.chr22.GT.FORMAT"
genotypes = fread(genotype_filename, header = T) # takes 4 minutes 

# subset to chr22: 
# genotypes_bak = genotypes
# genotypes = genotypes[CHROM = 'chr22',]

# remove missing genotypes and multi-allelic loci: 
genotypes = removeMissingGenotypes(genotypes)
genotypes = removeMultiAllelicLoci(genotypes)
dosages = genotypeToDosage(genotypes)

# calculate IBD2* and informative snp ratio for every pair of samples: 
samples = colnames(dosages)
ibs2_het_ratio_table = matrix(0, nrow = length(samples), ncol = length(samples))
colnames(ibs2_het_ratio_table) = rownames(ibs2_het_ratio_table) = samples
informative_snps_table = matrix(0, nrow = length(samples), ncol = length(samples))
colnames(informative_snps_table) = rownames(informative_snps_table) = samples

for (sample1 in samples){
	for (sample2 in samples){
		if (sample1 == sample2) next
		dosage1 = dosages[,sample1]
		dosage2 = dosages[,sample2]
		stopifnot(length(dosage1)==length(dosage2))
		ibs2_het = sum((dosage1 == dosage2) & (dosage1 == 1))
		ibs0 = sum(abs(dosage1 - dosage2) == 2)
		ibs1 = sum(abs(dosage1 - dosage2) == 1)
		ibs2_het_ratio = ibs2_het/(ibs0 + ibs2_het)
		ibs2_het_ratio_table[sample1, sample2] = ibs2_het_ratio
		informative_snps = (ibs0 + ibs2_het)/(ibs0 + ibs1 + ibs2_het)
		informative_snps_table[sample1, sample2] = informative_snps
	}
}

ibs2_het_ratio_table = as.data.frame(ibs2_het_ratio_table)
ibs2_het_ratio_table$sample1 = rownames(ibs2_het_ratio_table)
ibs2_het_ratio_table_long = melt(ibs2_het_ratio_table, id.vars = 'sample1', variable.name = 'sample2', value.name = 'ibs2_het_ratio')
ibs2_het_ratio_table_long$comparison = with(ibs2_het_ratio_table_long, paste(sample1, sample2, sep="_"))

informative_snps_table = as.data.frame(informative_snps_table)
informative_snps_table$sample1 = rownames(informative_snps_table)
informative_snps_table_long = melt(informative_snps_table, id.vars = 'sample1', variable.name = 'sample2', value.name = 'informative_snp')
informative_snps_table_long$comparison = with(informative_snps_table_long, paste(sample1, sample2, sep="_"))


# keep only comparison within same ethnic group: 
ancestry_proportion_filename = '../processed_data/ancestry_proportions.tsv'
ancestry_proportion = fread(ancestry_proportion_filename, header =T)
ancestry_proportion$genomic_ancestry = ancestry_proportion$self_report
ancestry_proportion[ID == '1051601', genomic_ancestry:= 'EUR'] 
ancestry_proportion[ID == '1448', genomic_ancestry:= 'AMR'] 
ancestry_proportion[ID == '150328', genomic_ancestry:= 'EUR'] 
ancestry_proportion[ID == '1848', genomic_ancestry:= 'AMR'] 
ancestry_proportion[ID == '1858', genomic_ancestry:= 'AMR'] 
ancestry_proportion[ID == '2102', genomic_ancestry:= 'AMR'] 
ancestry_proportion[ID == '24156', genomic_ancestry:= 'AFR'] 
ancestry_proportion[ID == '24635', genomic_ancestry:= 'AFR'] 
ancestry_proportion[ID == '289727', genomic_ancestry:= 'EUR'] 
ancestry_proportion[ID == '3100203', genomic_ancestry:= 'EUR'] 
ancestry_proportion[ID == '313605', genomic_ancestry:= 'AMR'] # only small AMR ancestry, but then this individual could only be AMR
ancestry_proportion[ID == '317155', genomic_ancestry:= 'AMR'] # only small AMR ancestry, but then this individual could only be AMR
ancestry_proportion[ID == '7103002', genomic_ancestry:= 'AFR'] 
ibs2_het_ratio_table_long$same_population = TRUE
ibs2_het_ratio_table_long$population = NA
informative_snps_table_long$same_population = TRUE
informative_snps_table_long$population = NA

for (i in 1:nrow(ibs2_het_ratio_table_long)){
	sample1 = ibs2_het_ratio_table_long$sample1[i]
	sample2 = ibs2_het_ratio_table_long$sample2[i]
	# print(sample1)
	# print(sample2)
	sample1_ancestry = ancestry_proportion[ID == sample1,genomic_ancestry]
	sample2_ancestry = ancestry_proportion[ID == sample2,genomic_ancestry]
	# print(sample1_ancestry)
	# print(sample2_ancestry)
	ibs2_het_ratio_table_long$population[i] = ifelse(sample1_ancestry==sample2_ancestry, sample1_ancestry,NA)
	informative_snps_table_long$population[i] = ifelse(sample1_ancestry==sample2_ancestry, sample1_ancestry,NA)

}
ibs2_het_ratio_table_long = as.data.table(ibs2_het_ratio_table_long)
informative_snps_table_long = as.data.table(informative_snps_table_long)
pdf()
plot(ibs2_het_ratio_table_long[!is.na(population), ibs2_het_ratio], informative_snps_table_long[!is.na(population), informative_snp])
dev.off()