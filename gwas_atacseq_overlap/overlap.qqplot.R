library(data.table)
library(dtplyr)
library(dplyr)
library(cowplot)
library(gap)
source('gwas_atacseq_overlap/utils.R')

# read GWAS data:
gwas=fread('/srv/persistent/bliu2/HCASMC_eQTL/data/gwas/cad.add.160614.website.txt')


# format GWAS data:
gwas=gwas%>%select(chrom=chr,pos=bp_hg19,rsid=markername,pval=p_dgc)
gwas=gwas%>%mutate(chromStart=pos,chromEnd=pos,chrom=paste0('chr',chrom))
gwas=gwas%>%select(chrom,chromStart,chromEnd,rsid,pval)


# read ATACseq peaks:
atac_hcasmc=readAtac('zcat ~/atacseq/2305/out/peak/idr/optimal_set/2305_ppr.IDR0.1.filt.narrowPeak.gz')
overlaps_hcasmc=overlap(atac_hcasmc,gwas)
qqunif(overlaps_hcasmc$pval)


# overlap GWAS and ATACseq:
atac_gastric=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E094-DNase.macs2.narrowPeak.gz')
overlaps_gastric=overlap(atac_gastric,gwas)
qq_gastric=qqunif(overlaps_gastric$pval,plot.it=F)
points(qq_gastric)


# 
atac_gastric=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E094-DNase.macs2.narrowPeak.gz')
overlaps_gastric=overlap(atac_gastric,gwas)
qq_gastric=qqunif(overlaps_gastric$pval,plot.it=F)
points(qq_gastric)


# 
atac_ovary=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E097-DNase.macs2.narrowPeak.gz')
overlaps_ovary=overlap(atac_ovary,gwas)
qq_ovary=qqunif(overlaps_ovary$pval,plot.it=F)
points(qq_ovary)


# 
atac_pancreas=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E098-DNase.macs2.narrowPeak.gz')
overlaps_pancreas=overlap(atac_pancreas,gwas)
qq_pancreas=qqunif(overlaps_pancreas$pval,plot.it=F)
points(qq_pancreas)


# 
atac_psoasMuscle=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E100-DNase.macs2.narrowPeak.gz')
overlaps_psoasMuscle=overlap(atac_psoasMuscle,gwas)
qq_psoasMuscle=qqunif(overlaps_psoasMuscle$pval,plot.it=F)
points(qq_psoasMuscle,col='red')


# 
atac_smallIntestine=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E109-DNase.macs2.narrowPeak.gz')
overlaps_smallIntestine=overlap(atac_smallIntestine,gwas)
qq_smallIntestine=qqunif(overlaps_smallIntestine$pval,plot.it=F)
points(qq_smallIntestine,col='green')


# check the score on each peak: 
hist(atac_smallIntestine$signalValue)
hist(atac_hcasmc$signalValue)
ggplot(atac_smallIntestine,aes(pValue))+geom_histogram()+geom_histogram(data=atac_hcasmc,aes(pValue),fill='red',alpha=0.2)
ggplot(atac_smallIntestine,aes(qValue))+geom_histogram()+geom_histogram(data=atac_hcasmc,aes(qValue),fill='red',alpha=0.2)
min(atac_hcasmc$pValue)
min(atac_smallIntestine$pValue)
nrow(atac_smallIntestine[pValue>min(atac_hcasmc$pValue)])
ggplot(atac_smallIntestine[pValue>min(atac_hcasmc$pValue)],aes(qValue))+geom_histogram()+geom_histogram(data=atac_hcasmc,aes(qValue),fill='red',alpha=0.2)
qqplot(atac_smallIntestine$pValue,atac_hcasmc$pValue)
nrow(atac_smallIntestine)


# rejection sampling:
ecdf_hcasmc=ecdf(atac_hcasmc$pValue)
ecdf_smallIntestine=ecdf(atac_smallIntestine$pValue)
keep=logical(nrow(atac_smallIntestine))
for (i in 1:nrow(atac_smallIntestine)){
	pval_sample=atac_smallIntestine$pValue[i]
	dsample=diff(ecdf_smallIntestine(c(pval_sample-0.5,pval_sample+0.5)))
	dtarget=diff(ecdf_hcasmc(c(pval_sample-0.5,pval_sample+0.5)))
	u=runif(1)
	keep[i]=(u<dtarget/dsample)
}
atac_smallIntestine_bak=atac_smallIntestine
atac_smallIntestine=atac_smallIntestine_bak[keep,]


# re-make qqplot: 
atac_smallIntestine=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E109-DNase.macs2.narrowPeak.gz')
overlaps_smallIntestine=overlap(atac_smallIntestine,gwas)
qq_smallIntestine=qqunif(overlaps_smallIntestine$pval,plot.it=F)
points(qq_smallIntestine,col='green')


# take the strongest gwas variant for each peak: 
overlaps_hcasmc_maxPval=overlaps_hcasmc%>%group_by(chrom,chromStart,chromEnd)%>%summarize(pval_max=max(pval))
overlaps_hcasmc_maxPval

overlaps_smallIntestine_maxPval=overlaps_smallIntestine%>%group_by(chrom,chromStart,chromEnd)%>%summarize(pval_max=max(pval))
overlaps_smallIntestine_maxPval


# re-make qqplot:
qqunif(overlaps_hcasmc_maxPval$pval)
qq_smallIntestine=qqunif(overlaps_smallIntestine_maxPval$pval,plot.it=F)
points(qq_smallIntestine,col='green')


# look at top hits: 
min(overlaps_smallIntestine_maxPval$pval)
min(overlaps_hcasmc_maxPval$pval)
overlaps_smallIntestine_maxPval[which(overlaps_smallIntestine_maxPval$pval==min(overlaps_smallIntestine_maxPval$pval))]
overlaps_smallIntestine[chrom=='chr9'&chromStart==22098546&chromEnd==22099083,]
# makr qqplot of GWAS variants in ATACseq regions:

