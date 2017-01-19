library(data.table)
library(cowplot)
library(dplyr)
library(VennDiagram)
library(limma)

in_file=commandArgs(T)[1]
fig_file=commandArgs(T)[2]
# in_file='../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_variants.allFeatures.bed'
# fig_file='../figures/mpra/eCAVIAR/qc_eCAVIAR_variants.pdf'


data=fread(in_file)
clpp_cut=0.01

# Plot effect size (pi):
data[,CLPP:=ifelse(clpp>clpp_cut,'CLPP>0.01','CLPP<0.01')]
p1=ggplot(data,aes(x=pi))+facet_grid(CLPP~.,scales='free_y')+geom_histogram(position='dodge',binwidth=0.01)+xlab('Estimated allelic ratio')


# Plot CADD score:
p2=ggplot(data,aes(x=caddPhred,color=CLPP))+geom_density(position='dodge')


# Plot ATACseq signal value: 
data[,ATACseqFill:=ifelse(is.na(ATACseq),0,ATACseq)]
p3=ggplot(data,aes(x=ATACseqFill))+facet_grid(CLPP~.,scales='free_y')+geom_histogram(position='dodge')


# Plot relationship between GWAS and eQTL posterior: 
p4=ggplot(data,aes(gwasPost,eqtlPost,color=CLPP))+geom_point(alpha=0.5,size=3)+xlab('GWAS posterior')+ylab('eQTL posterior')


# Plot correlation between prior and posterior genotype (test SNP): 
p5=ggplot(data,aes(rsq_rsnp,color=CLPP))+geom_density()+xlab('Correlation between prior and post. genotype (r2)')


# Plot PhastCons scores:
phastCons=melt(data,measure.var=c('priPhCons','mamPhCons','verPhCons'),value.name='phastCons',variable.name='type')
p6=ggplot(phastCons,aes(phastCons))+facet_grid(CLPP~type,scales='free_y')+geom_histogram(binwidth=0.01)


# Plot PhyloP scores: 
phyloP=melt(data,measure.var=c('priPhyloP','mamPhyloP','verPhyloP'),value.name='phyloP',variable.name='type')
p7=ggplot(phyloP,aes(phyloP))+facet_grid(CLPP~type,scales='free')+geom_histogram()


# Plot GERP scores: 
gerp=melt(data,measure.var=c('GerpN','GerpS'),value.name='Gerp',variable.name='type')
p8=ggplot(gerp,aes(Gerp))+facet_grid(CLPP~type,scale='free')+geom_histogram()


# Plot number of TFBS sites:
p9=ggplot(data,aes(n_tfbs))+facet_grid(CLPP~.,scale='free')+geom_histogram()+xlab('Number of TFBS')


# Plot CpG islands: 
p10=ggplot(data,aes(CpG))+facet_grid(CLPP~.,scale='free')+geom_histogram()


# Plot a Venn diagram for GWAS and eQTL putative variant sets.
vennCounts=data[,.(chr,end,gwasSet,eqtlSet)]%>%unique()%>%select(gwasSet,eqtlSet)%>%vennCounts()
colnames(vennCounts)=c('GWAS set','eQTL set','Counts')
vennDiagram(vennCounts)
p11=recordPlot()


# Plot effect size (pi):
data[,Set:=ifelse(gwasSet==1&eqtlSet==1,'Intersection','GWAS only')]
data[,Set:=ifelse(gwasSet==0&eqtlSet==1,'eQTL only',Set)]
p12=ggplot(data,aes(x=pi))+facet_grid(Set~.,scales='free_y')+geom_histogram(position='dodge',binwidth=0.001)+xlab('Estimated allelic ratio')


# Plot CADD score:
p13=ggplot(data,aes(x=caddPhred,color=Set))+geom_density(position='dodge')


# Plot ATACseq signal value: 
p14=ggplot(data,aes(x=ATACseqFill))+facet_grid(Set~.,scales='free_y')+geom_histogram(position='dodge')


# Plot relationship between GWAS and eQTL posterior: 
p15=ggplot(data,aes(gwasPost,eqtlPost,color=Set))+geom_point(alpha=0.5,size=3)+xlab('GWAS posterior')+ylab('eQTL posterior')


# Plot correlation between prior and posterior genotype (test SNP): 
p16=ggplot(data,aes(rsq_rsnp,color=Set))+geom_density()+xlab('Correlation between prior and post. genotype (r2)')+annotate(geom='text',x=0.5,y=50,label=paste0(data%>%filter(Set=='Intersection',rsq_rsnp<0.8)%>%nrow(),',',data%>%filter(Set=='GWAS only',rsq_rsnp<0.8)%>%nrow(),',',data%>%filter(Set=='eQTL only',rsq_rsnp<0.8)%>%nrow(),' variants in intersection, \nGWAS only and eQTL only with r2 < 0.8'))


# Plot PhastCons scores:
phastCons=melt(data,measure.var=c('priPhCons','mamPhCons','verPhCons'),value.name='phastCons',variable.name='type')
p17=ggplot(phastCons,aes(phastCons))+facet_grid(Set~type,scales='free_y')+geom_histogram(binwidth=0.01)


# Plot PhyloP scores: 
phyloP=melt(data,measure.var=c('priPhyloP','mamPhyloP','verPhyloP'),value.name='phyloP',variable.name='type')
p18=ggplot(phyloP,aes(phyloP))+facet_grid(Set~type,scales='free')+geom_histogram()


# Plot GERP scores: 
gerp=melt(data,measure.var=c('GerpN','GerpS'),value.name='Gerp',variable.name='type')
p19=ggplot(gerp,aes(Gerp))+facet_grid(Set~type,scale='free')+geom_histogram()


# Plot number of TFBS sites:
p20=ggplot(data,aes(n_tfbs))+facet_grid(Set~.,scale='free')+geom_histogram()+xlab('Number of TFBS')


# Plot CpG islands: 
p21=ggplot(data,aes(CpG))+facet_grid(Set~.,scale='free')+geom_histogram()



# Save plots: 
pdf(fig_file)
for (i in 1:21) {p=get(paste0('p',i)); print(p)}
dev.off()