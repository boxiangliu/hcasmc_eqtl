library(data.table)

# Variables:
gwas_fn='../data/gwas/CARDIoGRAMplusC4D/cad.add.160614.website.txt'


# Read GWAS data:
gwas=fread(gwas_fn)
gwas[,logp:=-log10(p_dgc)]
summary(gwas$logp)
gwas[which(gwas$p_dgc==min(gwas$p_dgc)),]
# Make QQ-plot:
qqnorm(gwas$logp)