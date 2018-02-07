library(dtplyr)
library(data.table)
library(stringr)
library(cowplot)
library(ggrepel)
library(MASS)
library(dplyr)


# command line inputs: 
args=commandArgs(T)
pmplot=args[1]
pval_vs_size=args[2]
# pmplot='../figures/tarid/tcf21_rs2329429.pmplot.pdf'
# pval_vs_size='../figures/tarid/tcf21_rs2329429.pvalue_vs_sample_size.pdf'


# read metasoft form stdin: 
metasoft=read.table(file("stdin"), header=F, fill=T, sep="\t")%>%as.data.table()
metasoft[,V107:=NULL]

# read tissue name (study name):
study_name=unlist(fread('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/Metasoft_tissue_order.alphabetical.txt',header=F))


# set metasoft column names:
col_names=c(paste('pvalue',study_name,sep='_'),paste('mvalue',study_name,sep='_'))
col_names=c(paste0('V',seq(1,16)),col_names)
stopifnot(ncol(metasoft)==length(col_names))
setnames(metasoft,col_names)


# subset to mvalues: 
mvalue=metasoft%>%dplyr::select(contains('mvalue'))%>%as.data.frame()
rownames(mvalue)=metasoft$RSID
colnames(mvalue)=str_replace(colnames(mvalue),'mvalue_','')


# subset to pvalues: 
pvalue=metasoft%>%dplyr::select(contains('pvalue'))%>%as.data.frame()
rownames(pvalue)=metasoft$RSID
colnames(pvalue)=str_replace(colnames(pvalue),'pvalue_','')


# read tissue names: 
tissue_names_file='/srv/persistent/bliu2/HCASMC_eQTL/scripts/160603/collapsed_tissue_names.3.txt'
tissue_names=fread(tissue_names_file,header=T)


# make PM plot: 
to_plot=data.frame(mvalue=unlist(mvalue),pvalue=-log10(unlist(pvalue)))
to_plot$tissue=rownames(to_plot)
to_plot_bak=to_plot
to_plot=merge(to_plot_bak,tissue_names,by.x='tissue',by.y='original')
set.seed(42)
p1=ggplot(to_plot,aes(mvalue,pvalue,label=tissue,color=collapsed))+geom_point(size=5)+theme_bw()+xlab('M-value')+ylab('-log10(P-value)')+geom_vline(xintercept=0.9,color='red',linetype='dashed')+scale_color_discrete(guide=F)+geom_text_repel(aes(mvalue,pvalue,label=ifelse(mvalue>0.9|tissue=='HCASMC',tissue,NA)),size=5,nudge_x=-0.2)
save_plot(pmplot,p1,base_width=6,base_height=6)



# make plot about sample size vs pvalue:
if (pval_vs_size!='none'){
	gtex_file='../processed_data/160530/gtex.v6p.egenes.summary.txt'
	gtex=read.table(gtex_file,sep='\t',header=T)%>%select(tissue=Tissue,size=Size)
	gtex=rbind(gtex,data.frame(tissue='HCASMC',size=52))
	to_plot=merge(to_plot,gtex,by='tissue')
	to_plot=to_plot%>%mutate(effect=ifelse(mvalue>0.9,'Yes','No'))
	set.seed(42)
	p2=ggplot(to_plot,aes(x=size,y=pvalue,color=effect))+geom_point(size=5)+stat_smooth(method='rlm',formula=y~x-1)+geom_text_repel(aes(x=size,y=pvalue,label=ifelse(effect=='Yes',tissue,NA)),show.legend=F)+theme_bw()+xlab('Size')+ylab('-log10(P-value)')
	save_plot(pval_vs_size,p2,base_height=6,base_width=6)
}
