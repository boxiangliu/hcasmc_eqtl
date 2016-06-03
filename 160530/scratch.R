pval_file='../processed_data/160530/cis.txt'
genotype_file='../processed_data/160530/dosage.tsv'

# read input:
pval=fread(pval_file,header=T)
pval_sample=pval[sample(nrow(pval),1000000),]
pval_sample_sort=sort(unlist(pval_sample[,"p-value",with=F]))
pval_sample_rank=rank(unlist(pval_sample[,"p-value",with=F]),ties.method='min')
pval_sample[which(pval_sample_rank==1)]
pval_sample[which(pval_sample_rank==3)]
length(unique(pval_sample_rank))
genotypes=fread(genotype_file,header=T)

genotypes[id%in%c('chr7_73243786_T_G','chr7_73255165_C_T')]
# these two variants have exactly the same genotypes



