library(data.table)
library(dplyr)
library(dtplyr)
library(stringr)
library(foreach)
library(doMC)
registerDoMC(10)
library(cowplot)

in_dir_released_all_cell_type='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_released_all_cell_type/'
pval_fn='../processed_data/gwas_atacseq_overlap/gregor/overlap_enrichment/gregor_pval_released_all_cell_type.tsv'
n_overlap_fn='../processed_data/gwas_atacseq_overlap/gregor/overlap_enrichment/gregor_num_overlap_released_all_cell_type.tsv'
out_dir='../processed_data/plot_n_peaks_vs_pval/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

# count the number of peaks:
fn_list=list.files(in_dir_released_all_cell_type,full.name=TRUE)

base_coverage=foreach(i=1:length(fn_list),.combine='rbind')%dopar%{

	sample=basename(fn_list[i])%>%str_replace(.,'.merged.bed','')%>%str_replace_all(.,' ','_')%>%str_replace_all(.,'\'','_')
	message(sample)

	cmd=sprintf('sort -k1,1 -k2,2n "%s" | bedtools merge > %s/%s.bed',fn_list[i],out_dir,sample)
	system(cmd)

	n_base=fread(sprintf('%s/%s.bed',out_dir,sample),
		col.names=c('chr','start','end'))%>%
	mutate(length=end-start)%>%
	summarize(n_base=sum(length))%>%
	unlist()

	data.table(sample,n_base)
}

# Read GREGOR p-value file:
gregor=fread(pval_fn)%>%mutate(tissue=str_replace_all(tissue,' ','_')%>%str_replace_all('\'',''))


# Merge base coverage and GREGOR p-value:
merged=merge(gregor,base_coverage,by.x='tissue',by.y='sample')


# Read the number of overlaps:
n_overlap=fread(n_overlap_fn)%>%mutate(sample=str_replace_all(sample,' ','_')%>%str_replace_all('\'',''))


# Merge with number of overlap: 
merged=merge(merged,n_overlap,by.x='tissue',by.y='sample')


# Plot the p-value vs number of peaks:
ggplot(merged[n_overlap>=2],aes(-log10(pval),n_base))+geom_point()+stat_smooth(method='lm')
ggplot(merged[n_overlap>=2],aes(n_base,-log10(pval)))+geom_point()+stat_smooth(method='lm')
nrow(merged[n_overlap>=5])
ggplot(merged[n_overlap>=2],aes(n_base,-log10(pval)))+geom_point()+stat_smooth(method='lm')+scale_x_log10()
fit=lm(pval~n_base,data=merged[n_overlap>=2])
summary(fit)