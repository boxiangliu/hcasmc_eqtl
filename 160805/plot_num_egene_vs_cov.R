#!/usr/bin/env Rscript
# bosh liu
# durga
# plot the number of egenes vs number of genotype PCs and number of PEER factors


# library
library(stringr)
library(gtools)
library(cowplot)


# command args:
figure_path='../figures/160805/num_egene_vs_cov.pdf'


# list input files:
input_files=list.files(path='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/find_optimum_num_PEER_factors',pattern='padj.txt',full.names=T)
input_files=mixedsort(input_files)


# read input: 
fastqtl=data.table()
for (input_file in input_files){
	input=fread(input_file,header=T)
	num_pc=str_split_fixed(basename(input_file),pattern="\\.",n=6)[,3] %>% str_replace('pc','') 
	num_peer=str_split_fixed(basename(input_file),pattern="\\.",n=6)[,4] %>% str_replace('peer','') %>% as.integer()
	input$num_pc=num_pc
	input$num_peer=num_peer
	fastqtl=rbind(fastqtl,input)
}


# calculate the number of egenes: 
num_egenes = fastqtl %>% group_by(num_pc,num_peer) %>% summarize(num_eqtls=sum(qval<0.05)) 
head(num_egenes)


# plot the number of egenes: 
p=ggplot(num_egenes,aes(x=as.factor(num_peer),y=num_eqtls,group=num_pc,color=num_pc))+geom_point()+geom_line()+theme(axis.text.x=element_text(angle=90,vjust=0.5))+xlab('Number of PEER factors')+ylab('Number of eGenes')
save_plot(figure_path,p)
