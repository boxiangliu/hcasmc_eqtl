library(data.table)
library(dtplyr)
library(dplyr)
library(cowplot)
source('gwas_atacseq_overlap/utils.R')


# 
gwas=fread('../processed_data/160615/compare_hcasmc_and_gtex/GWAS.txt')
setnames(gwas,c('chr','pos','rsid','id','type'))
gwas=gwas%>%mutate(chromStart=pos,chromEnd=pos)%>%select(chrom=chr,chromStart,chromEnd,rsid,id)
gwas[,chrom:=paste0('chr',chrom)]


# 
atac_hcasmc=readAtac('zcat ~/atacseq/2305/out/peak/idr/optimal_set/2305_ppr.IDR0.1.filt.narrowPeak.gz')
overlaps_hcasmc=overlap(atac_hcasmc,gwas)


# 
atac_gastric=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E094-DNase.macs2.narrowPeak.gz')
overlaps_gastric=overlap(atac_gastric,gwas)


# 
atac_ovary=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E097-DNase.macs2.narrowPeak.gz')
overlaps_ovary=overlap(atac_ovary,gwas)


# 
atac_pancreas=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E098-DNase.macs2.narrowPeak.gz')
overlaps_pancreas=overlap(atac_pancreas,gwas)


# 
atac_psoasMuscle=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E100-DNase.macs2.narrowPeak.gz')
overlaps_psoasMuscle=overlap(atac_psoasMuscle,gwas)


# 
atac_smallIntestine=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E109-DNase.macs2.narrowPeak.gz')
overlaps_smallIntestine=overlap(atac_smallIntestine,gwas)


# 
tissues=c('hcasmc','gastric','ovary','pancreas','psoasMuscle','smallIntestine')
overlap_summary=data.frame()
for (this_tissue in tissues){
	atac=get(paste0('atac_',this_tissue))
	overlaps=get(paste0('overlaps_',this_tissue))


	# Number of GWAS SNPs in each peak:
	assign(paste0('overlaps_peakSummary_',this_tissue),peakSummary(overlaps))


	# Number of unique peaks overlapping GWAS SNPs:
	n_unique_peaks=nrow(get(paste0('overlaps_peakSummary_',this_tissue)))


	tmp_df=data.frame(n_peaks=nrow(atac),n_overlaps=nrow(overlaps),n_unique_peaks=n_unique_peaks,tissue=this_tissue)
	overlap_summary=rbind(overlap_summary,tmp_df)
}
overlap_summary=as.data.table(overlap_summary)


#
overlap_summary[,fraction:=n_overlaps/n_peaks]


# 
p1=ggplot(overlap_summary,aes(x=reorder(tissue,-n_overlaps,FUN=mean),y=n_overlaps))+geom_bar(stat='identity')+xlab('Tissue')+ylab('Number of GWAS-DHS overlaps')
p2=ggplot(overlap_summary,aes(x=reorder(tissue,-fraction,FUN=mean),y=fraction))+geom_bar(stat='identity')+xlab('Tissue')+ylab('Fraction of GWAS-DHS overlaps')
p3=ggplot(overlap_summary,aes(x=reorder(tissue,-n_unique_peaks,FUN=mean),y=n_unique_peaks))+geom_bar(stat='identity')+xlab('Tissue')+ylab('Number of unique GWAS-DHS overlaps')
save_plot('../figures/gwas_atacseq_overlap/n_overlaps.pdf',p1,base_width=8,base_height=8)
save_plot('../figures/gwas_atacseq_overlap/fraction_overlaps.pdf',p2,base_width=8,base_height=8)
save_plot('../figures/gwas_atacseq_overlap/n_unique_overlaps.pdf',p3,base_width=8,base_height=8)