library(data.table)
library(stringr)
library(cowplot)


# Variable: 
fn='../processed_data/hcasmc_specific_open_chromatin/peak_specificity/HCASMC.bed'


# Read peak specificity data: 
peak=fread(fn,header=F,select=c(1:4,8,9),col.names=c('chr','start','end','signalValue','id','psi'))


# Parse the id column:
peak[,num_tissue:=as.integer(str_split_fixed(id,'_',4)[,4])]


# Assign subpeak within each peak a distance percentage: 
peak[,dist:=seq(1,.N)/(.N+1),by=c('chr','start','end')]


# Plot specificty score against the distance percentage: 
p1=ggplot(peak[sample(10000),],aes(dist,psi))+geom_point(alpha=0.01,size=5)+stat_smooth()+xlab('Position within peak')+ylab('Peak specificity score')


# Plot number of tissues sharing a subpeak against the distance percentage: 
p2=ggplot(peak[sample(10000),],aes(dist,num_tissue))+geom_point(alpha=0.01,size=5)+stat_smooth()+xlab('Position within peak')+ylab('Peak specificity score')


# Save plots: 
pdf('../figures/hcasmc_specific_open_chromatin/position_vs_specificity_and_num_tissues.pdf')
print(p1);print(p2)
dev.off()


# Number of peaks:
nrow(unique(peak[,.(chr,start,end)]))
# [1] 92681
# The original file in ../processed_data/hcasmc_specific_open_chromatin/\
# encode_plus_hcasmc_filt/HCASMC.bed has 92773 peaks. 


# Take the minimum of specificity score for each peak: 
peak[,min_psi:=min(psi),by=c('chr','start','end')]
peak[,dist:=NULL]
peak[,id:=NULL]
peak2=unique(peak[psi==min_psi,])
peak2[,min_psi:=NULL]


# Plot the distribution of specificity score:
pdf('../figures/hcasmc_specific_open_chromatin/hcasmc_specificity_distribution.pdf')
hist(peak2$psi,breaks=100,xlab='Specificity index',main='HCASMC')
dev.off()