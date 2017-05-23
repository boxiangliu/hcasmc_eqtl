#!/usr/bin/env Rscript
# bosh liu
# 2016/04/19
# durga 
# detect DNA sample contamination

args = commandArgs(T)
genotype_filename = args[1]
source('utils.R')

# genotype_filename = '../data/joint/recalibrated_variants.GT.FORMAT'
genotypes = fread(genotype_filename, header = T)

# subset to chr1, chr11, and chr22: 
genotypes_bak = genotypes
genotypes = genotypes[CHROM == 'chr22',]


# remove missing genotypes and multi-allelic loci: 
genotypes = removeMissingGenotypes(genotypes)
genotypes = removeMultiAllelicLoci(genotypes)
dosages = genotypeToDosage(genotypes)

het_hom_ratios = c()
for (sample in colnames(dosages)){
	sample_dosage = dosages[,sample]
	num_het = sum(sample_dosage == 1)
	num_nonref_hom = sum(sample_dosage == 2)
	het_hom_ratio = num_het/num_nonref_hom
	het_hom_ratios = c(het_hom_ratios, het_hom_ratio)
}

het_hom_ratios = data.frame(sample = colnames(dosages), ratio = het_hom_ratios)

# figure_path = '../figures/dna_contamination.pdf'
figure_path = args[2]
pdf(figure_path, width = 9, height = 6)
plot(het_hom_ratios, las =2, main = 'het/hom ratio')
text(het_hom_ratios, labels = het_hom_ratios$sample,srt = -45, adj = 0)
dev.off()
