source('utils.R')
args = commandArgs(T)
freq_filename = args[1]
strand_filename = args[2]
figure_filename = args[3]

# freq_filename = '../processed_data/imputation/recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.chr20.frq'
# strand_filename = '../processed_data/imputation/alignments.snp.strand'
# figure_filename = '../figures/plot_percentage_in_1000G.pdf'


freq = fread(freq_filename)
strand = fread(strand_filename)
setnames(strand, 'V1', 'TYPE')
setnames(strand, 'V4', 'SNP')

freq[, MAF := ifelse(MAF < 0.5, MAF, 1-MAF)] # correct MAF
freq[, BIN := cut(MAF, c(0,0.01,0.05,0.1,0.2,0.3,0.4,0.5))]
temp = merge(freq,strand[,.(TYPE,SNP)],by='SNP',all.x=T)
merged = temp

merged[,N_MISSING := sum(TYPE=='Missing',na.rm=T), by='BIN']
merged[,N_STRAND := sum(TYPE=='Strand',na.rm=T), by='BIN']
merged[,N := .N, by='BIN']
merged[,PCT_MISSING := N_MISSING/N, by='BIN']
merged[,PCT_STRAND := N_STRAND/N, by='BIN']

load('plot_percentage_in_1000G.RData')
p1 = ggplot(merged,aes(x=BIN)) + geom_point(aes(y=PCT_MISSING,color = 'PCT_MISSING')) + geom_point(aes(y=PCT_STRAND,color='PCT_STRAND')) + xlab('Frequency') + ylab('Percentage') + theme(axis.text.x = element_text(angle = 90, hjust=1)) + scale_y_continuous(breaks = seq(0,0.8,0.05))
save_plot(figure_filename, p1, base_width = 6)

