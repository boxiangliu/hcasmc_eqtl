# calculate the percentage of heterozygous sites are matched by its genome.
# the het sites can be from RNAseq, ATACseq, etc.
# 
# useage: 
# -sample	file list to sample in vcf format
# -reference	file to sample against in vcf format 

# setup:
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

# functions: 
cleanup <- function(dt) {
	dt = dt %>% select(1, 3, 4) %>% unique() 
	setnames(dt, c('chr','pos','ref_alt'))
	dt %>% group_by(chr, pos)
}

calculate_percentage_in_reference <- function(sample, reference){
	reference = reference %>% mutate(in_reference = 1)
	sample_reference = left_join(sample, reference, by = c('chr', 'pos', 'ref_alt')) 
	sample_reference[is.na(sample_reference)] = 0 
	percent_in_reference = mean(sample_reference$in_reference)
	
	# return;  
	percent_in_reference
}

# read command line args: 
options(echo = TRUE)
args = commandArgs(trailingOnly = TRUE)
sample_filename = args[1]; # if (str_detect(sample_filename, '.bed')) warning('input 1 should contain a list of bed files, one per line.')
reference_filename = args[2]
sample_filename = 'match_genomes.rnaseq.sample_list.txt'
reference_filename = '../data/joint/recalibrated_variants.GT.FORMAT.chr22'

# read wgs genotype:
reference = fread(reference_filename, header = T) 

# read rna genotypes: 
sample_filename = '../data/rnaseq/variants/chr22/1020301.var.flt.vcf'
sample = fread(sample_filename)
sample_filename_ls  = fread(sample_filename, header = FALSE, data.table = FALSE)
sample_ls = apply(sample_filename_ls, 1, function(x) cleanup(fread(x))) 
reference = cleanup(fread(input = reference_filename))

# calculate percent het sites from sample in reference: 
num_file = length(sample_ls)
percent_in_reference = data.frame(filename = rep(0, num_file), percent = rep(0, num_file))
for (i in seq(num_file)){
	sample = sample_ls[[i]]
	filename = sample_filename_ls[i,]
	percent_in_reference$percent[i] = calculate_percentage_in_reference(sample, reference)
	percent_in_reference$filename[i] = basename(filename)
}

# plot percent_in_reference:
pdf(file = '../figures/match_genome_percent_in_reference.pdf')
ggplot(percent_in_reference, aes(x = filename, y = percent)) + geom_bar(stat = 'identity') + ggtitle('percent het site match') + xlab('sample') + coord_flip()
dev.off()