#!/usr/bin/env Rscript
# bosh liu 
# 2016/04/06
# for each heterozygous sites in WGS, calculate the percentage of sites with 30-70% reference reads in the RNA variant count file
library(gplots)
source('utils.R')

args = commandArgs(T)
genotype_filename = args[1]
count_dir = args[2]
count_filenames = args[3]
figure_filename = args[4]
table_filename = args[5]
# genotype_filename = '../data/joint/recalibrated_variants.GT.FORMAT'
# count_dir = '../processed_data/rna_wgs_match/variant_count/'
# count_filenames = 'rna_wgs_match.R.sample_list.txt'
# figure_filename = '../figures/rna_wgs_match.pdf'

# read genotypes: 
genotypes = fread(genotype_filename, header = T)

# subset to chr1, chr11, and chr22: 
temp = genotypes[CHROM == 'chr1' | CHROM == 'chr11' | CHROM == 'chr22',]
genotypes = temp

# remove missing genotypes and multi-allelic loci: 
genotypes = removeMissingGenotypes(genotypes)
genotypes = removeMultiAllelicLoci(genotypes)
dosages = genotypeToDosage(genotypes)

# read count counts: 
count_files = scan(count_filenames, 'character')
counts = list()
for (count_file in count_files){
	counts[[length(counts) + 1]] = fread(paste0(count_dir, count_file))
	names(counts)[length(counts)] = dirname(count_file)
}


# for each individual in WGS, compare his/her genotype with all count files: 
stopifnot(nrow(dosages) == nrow(genotypes))
wgs_samples = colnames(dosages)
rna_samples = names(counts)
match_matrix = matrix(0, ncol = length(wgs_samples), nrow = length(rna_samples))
min_ref_ratio = 0.3
max_ref_ratio = 0.7
colnames(match_matrix) = wgs_samples
rownames(match_matrix) = rna_samples
for (wgs_sample in wgs_samples){
	message(paste("WGS:", wgs_sample))
	het_idx = which(dosages[,wgs_sample] == 1)
	het_loci = rownames(dosages)[het_idx]
	het_loci = as.data.table(data.frame(ID = het_loci))
	for (rna_sample in rna_samples){
		message(paste("RNAseq:",rna_sample))
		count = counts[[rna_sample]]
		count[,ID := paste(CHR, bp, sep = '_')]
		count_at_wgs_het = merge(het_loci, count, by = 'ID')
		percent_ref = count_at_wgs_het[,'%_REF',with=F]
		rna_het = which(percent_ref >= min_ref_ratio & percent_ref <= max_ref_ratio)
		percent_het = length(rna_het)/nrow(count_at_wgs_het)
		match_matrix[rna_sample, wgs_sample] = percent_het
	}
}

# make heatmap of match_matrix:  
pdf(figure_filename)
heatmap.2(match_matrix, trace = 'none', dendrogram='none', Rowv=F, Colv=F)
for (rna_sample in rna_samples){
	plot(match_matrix[rna_sample,], main = rna_sample, xlim=c(1,75), ylab = '% match', xlab = 'Sample index')
	text(match_matrix[rna_sample,], labels = wgs_samples, srt = -45, adj = -0.1)
}
for (wgs_sample in wgs_samples){
	plot(match_matrix[,wgs_sample], main = wgs_sample, xlim=c(1,84), ylab = '% match', xlab = 'Sample index')
	text(match_matrix[,wgs_sample], labels = rna_samples, srt = -45, adj = -0.1)
}
dev.off()

# make rna wgs hash table: 
best_matches = data.frame(rna = character(), wgs = character())
for (rna_sample in rna_samples){
	wgs_sample = names(which.max(match_matrix[rna_sample,]))
	best_matches=rbind(best_matches, data.frame(rna = rna_sample, dna=wgs_sample))
}
write.table(best_matches, file = table_filename, sep = '\t', row.names=FALSE, col.names=TRUE, quote=FALSE)