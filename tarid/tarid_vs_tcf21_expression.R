library(data.table)
library(cowplot)
library(dplyr)
library(dtplyr)


rpkm=fread('../processed_data/tarid/tarid_tcf21_rpkm.txt',header=T,check.name=F)
hash=fread('../data/gtex/SAMPID_SMTSD.sorted.filtered.with_hcasmc.txt',header=F)
setnames(hash,c('sample_ID','tissue'))
rpkm[,Name:=NULL]
rpkm[,Description:=NULL]

rpkm_bak=rpkm
rpkm=t(rpkm_bak)
rpkm=data.table(rpkm,keep.rownames=T)
setnames(rpkm,c('sample_ID','TARID','TCF21'))
rpkm=merge(rpkm,hash,by='sample_ID')
rpkm_median=rpkm%>%group_by(tissue)%>%summarize(TARID=median(TARID),TCF21=median(TCF21))

p1=ggplot(rpkm,aes(x=TARID,y=TCF21,color=tissue,group=1))+geom_point()+theme(legend.position="none")+stat_smooth(formula=y~x,method='lm')+scale_x_log10()+scale_y_log10()
p2=ggplot(rpkm_median,aes(x=TARID,y=TCF21,label=tissue,group=1,color=ifelse(tissue=='HCASMC',1,0)))+geom_point()+theme(legend.position="none")+stat_smooth(formula=y~x,method='lm')+scale_x_log10()+scale_y_log10()+geom_text(angle=-45)
save_plot('../figures/tarid/tarid_vs_tcf21.pdf',p1)
save_plot('../figures/tarid/tarid_vs_tcf21.median.pdf',p2,base_width=8,base_height=8)