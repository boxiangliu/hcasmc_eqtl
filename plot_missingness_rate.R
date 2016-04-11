library(data.table)

args = commandArgs(T)
# inp_prefix = '../data/joint/recalibrated_variants.GRCh37.pass'
# figure_path = '../figures/missingness.pdf'
inp_prefix = args[1]
figure_path = args[2]

imiss_filename = paste0(inp_prefix, '.imiss')
lmiss_filename = paste0(inp_prefix, '.lmiss')
IMISS=fread(imiss_filename, header=T)
LMISS=fread(lmiss_filename, header=T)
# subset for autosome: 
LMISS = LMISS[CHR!=23 & CHR!=24,]
LMISS_sample = LMISS[sample(nrow(LMISS),10000)]
pdf(figure_path, width = 14, height = 7)
oldpar=par(mfrow=c(1,2))
plot( (1:dim(IMISS)[1])/(dim(IMISS)[1]-1), sort(1-IMISS$F_MISS), main="Individual call rate cumulative distribution", xlab="Quantile", ylab="Call Rate" ); grid()
plot( (1:dim(LMISS_sample)[1])/(dim(LMISS_sample)[1]-1), sort(1-LMISS_sample$F_MISS), main="SNP coverage cumulative distribution", xlab="Quantile", ylab="Coverage" ); grid()
par(oldpar)
dev.off()