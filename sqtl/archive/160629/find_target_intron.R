#!/usr/bin/env Rscript
# durga
# bosh liu
# find putative target introns of 7 GWAS loci


# load sQTL dataset: 
sqtl=fread('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160629/sqtl.nominal.allpairs.normal.1e5.padj.txt',header=T)


# find genes influenced by GWAS loci:
# CDKN2BAS rs1537373 chr9:22103341
sqtl[geno=='chr9_22103341_T_G'&pval<0.1,.(geno,pheno,pval)]
# nothing


# SMAD3 rs17293632 chr15_67442596_C_T
sqtl[geno=='chr15_67442596_C_T'&pval<0.1, .(geno,pheno,pval)]
#                  geno                             pheno      pval
# 1: chr15_67442596_C_T chr15:67418322:67419633:clu_10853 0.0018936
# 2: chr15_67442596_C_T chr15:67418322:67457233:clu_10853 0.0561846
# 3: chr15_67442596_C_T chr15:67457426:67457591:clu_10854 0.0535347
# 4: chr15_67442596_C_T chr15:67457722:67459117:clu_10854 0.0346901



# PDGFD rs2019090 chr11:103668962
sqtl[geno=='chr11_103668962_A_T'&pval<0.1, .(geno,pheno,pval)]
# nothing


# IL6R rs7549250 chr1_154404336_C_T:
sqtl[geno=='chr1_154404336_C_T'&pval<0.1, .(geno,pheno,pval)]
# nothing

# LMOD1 rs34091558 1:201886770 
sqtl[str_detect(geno,'chr1_201886770')&pval<0.1, .(geno,pheno,pval)]
sqtl[str_detect(geno,'chr1_201886770'), .(geno,pheno,pval)]

# nothing

# BMP1 rs73551707 chr8_22046423_C_T
sqtl[geno=='chr8_22046423_C_T'&pval<0.1, .(geno,pheno,pval)]
#                 geno                            pheno      pval
# 1: chr8_22046423_C_T chr8:22054933:22056535:clu_32054 0.0340183


# CCDC97 rs2241718 chr19:41829606
sqtl[geno=='chr19_41829606_G_A'&pval<0.1, .(geno,pheno,pval)]
#                  geno                             pheno      pval
# 1: chr19_41829606_G_A chr19:41744059:41745069:clu_13762 0.0423558
# 2: chr19_41829606_G_A chr19:41770597:41774128:clu_13765 0.0422596