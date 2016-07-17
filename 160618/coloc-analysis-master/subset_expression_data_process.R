library(tidyr)
library(dplyr)
library(data.table)

map <- fread("/srv/scratch/tshimko/analysis/ExpressionData/sample_tissue_mappings.txt")
rpkm <- fread("/srv/scratch/tshimko/analysis/ExpressionData/sample_rpkm.txt")

long_rpkm <- rpkm %>% gather(SAMPID, rpkm, 3:ncol(rpkm))

output <- left_join(long_rpkm, map)

write.table(output$SAMPID, "tmp.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
subject <- fread("cat tmp.txt | tr '-' '\t' | cut -f1,2 | tr '\t' '-'", header = FALSE)
system("rm tmp.txt")

output$subject <- subject$V1

output <- output %>% select(Name, SAMPID, subject, SMTSD, rpkm) %>%
  rename(gene = Name, sample = SAMPID, subject = subject, tissue = SMTSD, rpkm = rpkm)

individuals <- fread("/srv/scratch/tshimko/analysis/ExpressionData/subset_tissues_individuals.txt")

subset_expression <- inner_join(output, individuals)

median_expression <- subset_expression %>% group_by(gene, tissue) %>% summarise(median_rpkm = median(rpkm))

write.table(median_expression, file = "/srv/scratch/tshimko/analysis/ExpressionData/median_expression_subset.txt", row.names = FALSE, quote = FALSE, sep = "\t")

write.table(output, file = "/srv/scratch/tshimko/analysis/ExpressionData/mapped_rpkm.txt", row.names = FALSE, quote = FALSE, sep = "\t")
