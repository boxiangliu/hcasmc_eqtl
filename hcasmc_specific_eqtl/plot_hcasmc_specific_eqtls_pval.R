library(data.table)
library(dplyr)
library(dtplyr)
library(cowplot)

eqtl=fread('../processed_data/eqtl_and_atacseq/specificity.mean.sorted.10000.eql.txt')
setnames(eqtl,c('sid','pval','beta','se','tissue'))
pdf('../figures/hcasmc_specific_eqtl/qsiTop1e5.pval.pdf')
ggplot(eqtl[tissue=='HCASMC',],aes(-log10(pval)))+geom_histogram()
dev.off()



