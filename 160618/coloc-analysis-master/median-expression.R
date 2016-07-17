library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
name <- args[1]
file <- paste0("/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/rna-seq/",
               args[1])
data <- fread(file) %>% as.data.frame(.)
melted_data <- gather(data, Individual, Expression, -Gene)
output <- melted_data %>% group_by(Gene) %>% summarise(MedianExp = median(Expression)) %>%
  mutate(Tissue = str_split(name, "\\.")[[1]][1])

new_name <- paste0(str_split(name, "\\.")[[1]][1], "_Expression.txt")

new_file <- paste0("/srv/scratch/tshimko/ExpressionLevels/", new_name)

write.table(output, file = new_file, sep = "\t", row.names = FALSE, quote = FALSE)
