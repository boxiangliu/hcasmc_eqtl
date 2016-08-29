pval=fread('../processed_data/160824/gene_min_pval.txt')
writeLines(pval[pval$V2<1e-4&pval$V3<1e-4,V1])


