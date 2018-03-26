library('cowplot')
library('ggrepel')
library('openxlsx')
library('data.table')
in_file='../processed_data/genotype/genotype_pc/genotype_pcs.52samples.tsv'
sample_sheet_file = '../data/sample_info/sample_info.xlsx'
fig_dir = '../figures/genotype/genotype_pc/'

sample_sheet=as.data.table(readWorkbook(sample_sheet_file,sheet=5))

# Make final figure:
genotype_pcs=read.table(in_file,check.names=FALSE)
data = as.data.frame(t(genotype_pcs[1:4,]))
data$ID = rownames(data)
data = merge(data,sample_sheet[,list(ID=RNA,Reported_Ethnicity)],by='ID')

p = ggplot(data,aes(C1,C2,label=ID,color=Reported_Ethnicity))+
	geom_point(size=3,alpha=0.5)+
	scale_color_discrete(name='Self-Reported Ancestry')+
	xlab('PC1')+
	ylab('PC2')

save_plot(sprintf('%s/pc_52samples_color_by_self_reported_ancestry.pdf',fig_dir),p,base_width=6.5)