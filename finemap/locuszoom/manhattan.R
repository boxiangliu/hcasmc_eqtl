library(data.table)
library(cowplot)
library(ggrepel)

smr_fn='../processed_data/finemap/smr/UKBB_HCASMC_smr_results.out.smr'
fm_fn='/users/mgloud/projects/brain_gwas/output/2017-11-13_08-51-41_hcasmc_full/manhattan/UKBB_GWAS1KG_EXOME_CAD_SOFT_META_PublicRelease_300517_txt_gz_hcasmc_eqtls_txt_gz.txt'
fig_dir='../figures/finemap/locuszoom/manhattan/'
if(!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}

smr=fread(smr_fn)[,list(chr=ProbeChr,gene_name=Gene,pos=Probe_bp,p=-log10(p_SMR),method='SMR')]
fm=fread(fm_fn)[,list(chr=chrom,pos,gene_name=gene,p=-clpp_score,method='eCAVIAR')]



data=rbind(smr,fm)
data[,color:=ifelse((chr%%2)&method=='SMR','color1',NA)]
data[,color:=ifelse((!chr%%2)&method=='SMR','color2',color)]
data[,color:=ifelse((chr%%2)&method=='eCAVIAR','color3',color)]
data[,color:=ifelse((!chr%%2)&method=='eCAVIAR','color4',color)]
data[,method:=factor(method,level=c('SMR','eCAVIAR'))]
data[,label:=ifelse(gene_name%in%c('FES','SMAD3','SIPA1','TCF21'),gene_name,'')]

dummy=data[,list(range=range(pos)),by='chr']
window=2e7
dummy[,pos:=ifelse(range==min(range),range-window,range+window),by=chr]
dummy[,c('p','color','label'):=list(0,NA,'')]



pdf(sprintf('%s/manhattan.pdf',fig_dir),height=4,width=8)
ggplot(data,aes(pos,p,color=color))+
	geom_point()+
	facet_grid(method~chr,scales='free',
		space='free_x',
		labeller=labeller(chr=function(x) {ifelse(x%in%c(1:19,21),x,'')}))+
	theme(axis.text.x=element_blank(),
		axis.title.x=element_blank(),
		axis.line.x=element_blank(),
		axis.ticks.x=element_blank(),
		panel.spacing=unit(0,'lines'),
		panel.border=element_blank(),
		strip.background=element_blank())+
	geom_blank(data=dummy)+
	scale_y_continuous(labels=function(x){abs(x)})+
	ylab(paste('CLPP                 -log10(P)'))+
	scale_color_manual(values=c(color1='grey',
		color2='black',color3='grey',color4='black'),guide=FALSE)
dev.off()


