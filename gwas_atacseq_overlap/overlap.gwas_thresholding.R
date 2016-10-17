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
atac_gastric=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E094-DNase.macs2.narrowPeak.gz')
atac_ovary=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E097-DNase.macs2.narrowPeak.gz')
atac_pancreas=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E098-DNase.macs2.narrowPeak.gz')
atac_psoasMuscle=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E100-DNase.macs2.narrowPeak.gz')
atac_smallIntestine=readAtac('zcat /mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/E109-DNase.macs2.narrowPeak.gz')



# overlap summary: 
thd=c(1e-3,1e-4,1e-5,1e-6,1e-7)
tissues=c('hcasmc','gastric','ovary','pancreas','psoasMuscle','smallIntestine')
overlap_summary=data.frame()
for (this_thd in thd){
	# thresholding gwas by pvalue: 
	gwas_thd=gwas[pval<this_thd,]


	for (this_tissue in tissues){
		# get atacseq data: 
		atac=get(paste0('atac_',this_tissue))


		# calculate overlaps: 
		overlaps=overlap(atac,gwas_thd)


		# Number of GWAS SNPs in each peak:
		assign(paste0('overlaps_peakSummary_',this_tissue),peakSummary(overlaps))


		# Number of unique peaks overlapping GWAS SNPs:
		n_unique_peaks=nrow(get(paste0('overlaps_peakSummary_',this_tissue)))


		# append to overlap summary: 
		tmp_df=data.frame(n_peaks=nrow(atac),n_overlaps=nrow(overlaps),n_unique_peaks=n_unique_peaks,tissue=this_tissue,threshold=this_thd)
		overlap_summary=rbind(overlap_summary,tmp_df)
	}
}


# calculating fractions:
overlap_summary=as.data.table(overlap_summary)
overlap_summary[,fraction_uniq:=n_unique_peaks/n_peaks]
overlap_summary[,fraction:=n_overlaps/n_peaks]


# calculating lfc:
overlap_summary[,rank:=rank(fraction_uniq),by='threshold']
overlap_summary=overlap_summary%>%arrange(threshold,desc(rank))
lfc=log10(overlap_summary[rank==max(rank),fraction_uniq]/overlap_summary[rank==max(rank)-1,fraction_uniq])
annotation=data.frame(overlap_summary[rank==max(rank),.(threshold,fraction_uniq,tissue)],lfc=prettyNum(lfc,digits=3))


# plotting: 
# p=ggplot(overlap_summary,aes(x=threshold,y=fraction_uniq,color=tissue))+geom_point()+scale_x_log10(breaks=10^(-seq(3,7)))+scale_y_log10(breaks=10^(-seq(2,4,by=0.3)))+theme_bw()+geom_text(data=annotation,aes(x=threshold,y=fraction_uniq,label=lfc),vjust=-1,color='black')+xlab('Threshold')+ylab('Fraction of unique overlaps')
p=ggplot(overlap_summary,aes(x=threshold,y=fraction_uniq,color=tissue))+geom_point()+scale_x_log10(breaks=10^(-seq(3,7)))+scale_y_log10(breaks=10^(-seq(2,4,by=0.3)))+theme_bw()+xlab('Threshold')+ylab('Fraction of unique overlaps')
save_plot('../figures/gwas_atacseq_overlap/gwas_atacseq_overlap.unique.pdf',p,base_height=4,base_width=6)

p1=ggplot(overlap_summary[threshold==1e-7,],aes(x=reorder(tissue,-fraction_uniq,FUN=mean),y=fraction_uniq,color=tissue))+geom_point()+scale_y_log10(breaks=10^(-seq(2,4,by=0.3)))+theme_bw()+xlab('Tissue')+ylab('Fraction of unique overlaps')+theme(axis.text.x=element_text(angle=45,hjust=1))
save_plot('../figures/gwas_atacseq_overlap/gwas_atacseq_overlap.unique.1e-7.pdf',p1,base_height=4,base_width=6)

