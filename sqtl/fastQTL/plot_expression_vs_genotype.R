#!/usr/bin/env Rscript
# bosh liu
# durga
# plot the expression vs genotype

# library:
library(cowplot)
library(data.table)


# command args:
genotype_file = '../processed_data/sqtl/fastQTL/plot_expression_vs_genotype/dosage.tsv'
expression_file = '../processed_data/sqtl/fastQTL/plot_expression_vs_genotype/expression.tsv'
fig_dir = '../figures/sqtl/fastQTL/plot_expression_vs_genotype/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

# read input:
expression=read.table(expression_file,header=T,row.names=1,check.names=FALSE)
genotype=read.table(genotype_file,header=T,row.names=1,check.names=FALSE)
setnames(expression,c('317155','2999','1508','1401','2105'),c('313605','289727','2109','CA1401','2102'))

data=data.frame(expression=unlist(expression))
data$genotype=as.character(unlist(genotype[,match(rownames(data),colnames(genotype))]))

p=ggplot(data,aes(genotype,expression))+
	geom_boxplot()+
	geom_jitter(width=0.2)

save_plot(sprintf('%s/19:44244370:44248924:clu_12328_rs4760.pdf',fig_dir),p)