library(data.table)
library(dtplyr)
library(dplyr)
library(cowplot)
library(gap)
library(stringr)
source('gwas_atacseq_overlap/utils.R')

# Variables:
in_dir='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_adult_filt/'
fig_dir='../figures/gwas_atacseq_overlap/gwas_thresholding/'
out_dir='../processed_data/gwas_atacseq_overlap/gwas_thresholding/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir)}
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# read GWAS data:
gwas=fread('../data/gwas/CARDIoGRAMplusC4D/cad.add.160614.website.txt')


# format GWAS data:
gwas=gwas%>%select(chr,pos=bp_hg19,rsid=markername,pval=p_dgc)
gwas=gwas%>%mutate(start=pos,end=pos,chr=paste0('chr',chr))


# create thresholded GWAS sets:
threshold=10^-seq(3,15,1)
gwas_ls=list()
for (i in 1:length(threshold)){
	t=threshold[i]
	gwas_ls[[i]]=gwas[pval<t,]
}


# overlap: 
fn=list.files(in_dir)
overlap_df=data.frame()
for (f in fn){
	tissue=str_replace(f,'.merged.bed','')
	print(sprintf('INFO - %s',tissue))
	x=fread(sprintf('%s/%s',in_dir,f))
	x=x[chr%in%paste0('chr',1:22)]
	setkey(x,chr,start,end)

	for (i in 1:length(threshold)){
		# thresholding gwas by pvalue: 
		t=threshold[i]
		gwas_thd=gwas_ls[[i]]
		setkey(gwas_thd,chr,start,end)

		# calculate overlaps:
		overlaps=foverlaps(x,gwas_thd,nomatch=0)
		tmp=nrow(unique(overlaps[,list(chr,pos,rsid,pval)]))
		overlap_df=rbind(overlap_df,data.frame(tissue,threshold=t,n_overlaps=tmp,n_peaks=nrow(x)))
	}
}

# Save:
setDT(overlap_df)
fwrite(overlap_df,sprintf('%s/num_overlap.tsv',out_dir),sep='\t')

# Plot:
pdf(sprintf('%s/num_overlap.pdf',fig_dir))
ggplot(overlap_df,aes(threshold,n_overlaps,color=ifelse(tissue=='HCASMC','red','grey')))+geom_point(alpha=0.5)+scale_color_discrete(guide='none')+scale_alpha_discrete(guide='none')+scale_x_log10()+scale_y_log10()+stat_smooth()
ggplot(overlap_df,aes(threshold,n_overlaps,color=ifelse(tissue=='HCASMC','red','grey'),group=tissue))+geom_point(alpha=0.5)+geom_line(alpha=0.5)+scale_color_discrete(guide='none')+scale_alpha_discrete(guide='none')+scale_x_log10()+scale_y_log10()
ggplot(overlap_df[tissue!='HCASMC',],aes(threshold,n_overlaps,group=as.character(threshold)))+geom_boxplot()+geom_jitter(alpha=0.5,height=0,width=0.1)+scale_x_log10()+ylim(-1,1100)+scale_y_log10()+geom_point(data=overlap_df[tissue=='HCASMC',],aes(threshold,n_overlaps),color='red',size=5,alpha=0.5)
dev.off()
