library(data.table)
library(stringr)
library(cowplot)

in_dir='../processed_data/verifyBamID/'
in_fn=list.files(in_dir,pattern='bestSM',full.names=TRUE)
fig_dir='../figures/rna_wgs_match/merge_bestSM_files/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

container=list()
for (i in 1:length(in_fn)){
	sample=str_extract(in_fn[i],pattern='(?<=/)([0-9]+)(?=\\.bestSM)')
	message('INFO - ', sample)
	temp=fread(in_fn[i])
	temp=data.table(temp,DNA=sample)
	container[[sample]]=temp
}

sample=Reduce(rbind,container)
stopifnot(all(sample[,`#SEQ_ID`==DNA]))

p=ggplot(sample,aes(x=as.character(`#SEQ_ID`),y=FREEMIX))+
	geom_point()+theme(axis.text.x=element_text(angle=90),axis.title.x=element_blank())+
	ylab('Contamination')

save_plot(sprintf('%s/freemix.pdf',fig_dir),p,base_width=7)