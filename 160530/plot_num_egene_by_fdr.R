#!/usr/bin/env Rscript 
# bosh liu
# plot the number and the proportion of eGenes by FDR:

# library:
library(cowplot)


# command args: 
args=commandArgs(T)
fastqtl_file=args[1]
fastqtl_file='../processed_data/160530/fastqtl.padj.txt'
figure_file=args[2]


# read input: 
fastqtl=fread(fastqtl_file,header=T)


# calculate the number of eGenes by FDR: 
fdrs=c(1e-3,0.01,0.05,0.1)
num_egene=c()
for (fdr in fdrs){
	num_egene=c(num_egene,sum(fastqtl$bh <= fdr))
}
to_plot=data.frame(fdr=fdrs,num_egene=num_egene)

# plot the number of eGenes by FDR: 
p=ggplot(to_plot,aes(x=as.factor(fdr),y=num_egene))+geom_bar(stat="identity")+xlab('FDR (BH)')+ylab('Number of eGenes')
p=p+geom_text(label=num_egene,nudge_y=15)
save_plot(figure_file,p)
