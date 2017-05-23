source('utils.R')
library(gplots)

# read genotypes: 
inp = "../data/joint/recalibrated_variants.GRCh37.biallelic.nomissing.GT.FORMAT"
genotypes = fread(inp, header = T)
stopifnot(ncol(genotypes) == 69)
genotypes = removeMissingGenotypes(genotypes)
detectMultiAllelicLoci(genotypes) # remove multi allelic loci
dosages = genotypeToDosage(genotypes) # convert genotypes to dosages
dim(dosages)

# make heatmap based on genotype correlation: 
correlation = cor(dosages)
figure1 = '../figures/genotype_correlation_heatmap.pdf'
pdf(figure1, width = 10, height = 10)
heatmap.2(correlation, trace = 'none')
dev.off()
correlation['2102','2105'] # 0.9855898
correlation['150328','59386145'] # 0.9998804
correlation['2999','289727'] # 0.9862771
correlation['1346','CA1346'] # 0.9859308
correlation['2109','CA1508'] # 0.9862455
correlation['317155','313605'] # 0.9870235

# compare 150328 and 59386145 (they should have the same genotype): 
genotype_150328 = dosages[,'150328']
genotype_59386145 = dosages[,'59386145']
diff = which(genotype_150328 != genotype_59386145)
length(diff)/length(genotype_150328) # 7.218316e-05

# output variants different between 150328 and 59386145 to bed file: 
diff_split = str_split_fixed(names(diff), "_", n = 2) 
diff_bed = data.frame(chrom = diff_split[,1], chromStart = as.integer(diff_split[,2])-1, chromEnd = as.integer(diff_split[,2]))
write.table(diff_bed, file = '../processed_data/variants_different_between_150328_and_59386145.bed', quote=FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

# load saved image: 
load('../r_environments/plot_genotype_correlation.RData')