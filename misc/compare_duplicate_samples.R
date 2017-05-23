source('utils.R')

args = commandArgs(T)
inp = args[1]
output1 = args[2]
output2 = args[3]


# read genotypes: 
genotypes = fread(inp)


# remove "GT" from column names:  
old = names(genotypes)
new = str_replace(old, '.GT', '')
setnames(genotypes, old, new)


# annotate indel or insertion deletion: 
genotypes[,TYPE := ifelse(nchar(REF)==1 & nchar(ALT)==1,'SNP','INDEL')]


# set missing variant ids:
ID = genotypes[,ID]
new_ID = genotypes[,paste(CHROM, POS, REF, ALT, sep="_")]
ID = ifelse(ID == '.', new_ID, ID)
genotypes$ID = ID


# duplicate samples: 
pairs = list()
pairs[[1]] = c('150328','59386145')
pairs[[2]] = c('2999','289727')
pairs[[3]] = c('1346','CA1346')
pairs[[4]] = c('2109','CA1508')
pairs[[5]] = c('317155','313605')
pairs[[6]] = c('2102','2105')


# subset to duplicate samples: 
all_pairs = unlist(pairs)
genotypes_dup = genotypes[,all_pairs,with=F]


# compare the difference between pairs: 
concordance = data.table()
for (pair in pairs){
	dup1 = genotypes_dup[,pair[1],with=F]
	dup2 = genotypes_dup[,pair[2],with=F]
	diff = (dup1 != dup2)
	colnames(diff) = NULL
	concordance = rbind(concordance, data.frame(CHROM = genotypes$CHROM, POS = genotypes$POS, diff = diff, pair = paste(pair[1], pair[2], sep = "_")))

} 
concordance = data.table::dcast(concordance, CHROM + POS ~ pair, value.var = 'diff')
concordance$num_diff = rowSums(concordance[,3:ncol(concordance),with=F])


# distribution of discordant loci:
tab = table(concordance$num_diff, genotypes$TYPE)
write.table(tab, file = output1, quote=F,sep='\t',row.names=T,col.names=T)



# output IDs for loci with more than 2 discordances: 
discordant_loci_ID = genotypes[concordance$num_diff > 2,ID]
write.table(discordant_loci_ID, file = output2, quote=F, row.names = F, col.names=F)