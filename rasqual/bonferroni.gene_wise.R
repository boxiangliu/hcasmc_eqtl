library(data.table)
library(qvalue)
pval=fread('../data/eQTL/rasqual/expressedGenes.pval.txt')
pval[,bonf:=p.adjust(pval,method='bonferroni'),by='fid']
pval[,min_bonf:=min(bonf),by='fid']
pval2=unique(pval[,.(fid,min_bonf)])
qval=qvalue(pval2$min_bonf)$qvalue
sum(qval<0.05) # 1438