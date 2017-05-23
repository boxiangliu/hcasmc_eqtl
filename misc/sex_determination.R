#!/usr/bin/env Rscript
# bosh liu
# 2016/04/16
# durga
# determines sex of individiuals

library(XLConnect)
source('utils.R')

args=commandArgs(T)
genotype_filename = args[1]
# genotype_filename = '../data/joint/recalibrated_variants.GT.FORMAT'
genotypes = fread(genotype_filename,header=T)

# subset to chromosome Y: 
chrY = genotypes[CHROM=='chrY',]

# calculate missingness rate on chrY:
missingY = (chrY[,-c('CHROM','POS'),with=F] == "./.")
missingY_pct = colMeans(missingY)

# plot missingness rate on chrY:
# figure_path = '../figures/sex.pdf'
figure_path = args[2]
pdf(figure_path)
plot(missingY_pct,ylab = "% missing on chrY", xlim = c(1,72), ylim = c(0,0.5), col = ifelse(missingY_pct > 0.3, 'red','blue'), main = 'chrY')
text(missingY_pct,names(missingY_pct),srt=45,adj= -0.1)
legend('right',legend = c('male','female'), col = c('blue','red'),pch = 1)


# subset to chromosome X: 
chrX = genotypes[CHROM=='chrX'] # 885778 rows 
chrX = removeMultiAllelicLoci(chrX) # 694245 rows 
chrX = removeMissingGenotypes(chrX) # no missing genotypes, because removeMultiAllelicLoci already removed them.
dosagesX = genotypeToDosage(chrX)

# calculate percentage of heterozygous sites on chrX:
heterozygousX = (dosagesX == 1)
heterozygousX_pct = colMeans(heterozygousX)

# plot the percentage of heterozygous sites on chrX:
plot(heterozygousX_pct, ylab = '% heterozygous on chrX', xlim = c(1,72), ylim = c(0, 0.35), col = ifelse(heterozygousX_pct > 0.05, 'red','blue'), main = 'chrX')
text(heterozygousX_pct,names(heterozygousX_pct),srt=45,adj= -0.1)
legend('topright',legend = c('male','female'), col = c('blue','red'),pch = 1)
dev.off()


# create data.frame for sex: 
stopifnot(all(names(missingY_pct) == names(heterozygousX_pct)))
sex = data.frame(sample= names(missingY_pct), chrY = ifelse(missingY_pct > 0.3, 'F','M'), chrX =  ifelse(heterozygousX_pct > 0.05, 'female','male'))

# add reported sex: 
sex = read.table(table_path,header=T,as.is=T)
sample_sheet_filename = '../processed_data/R21R33_Cell_line_summary.xlsx'
sample_sheet = readWorksheet(loadWorkbook(sample_sheet_filename),sheet=1)
merged = merge(sex, sample_sheet[,c('ID','Gender')], by.x='sample', by.y = 'ID')
setnames(merged, 'Gender','reported')
merged$concordant = (merged$chrY == merged$reported) & (merged$chrX == merged$reported)
merged
# save table: 
table_path = args[3]
table_path = '../processed_data/sex.tsv'
write.table(merged, file = table_path, quote = F, sep = '\t', row.names=F, col.names=T)
