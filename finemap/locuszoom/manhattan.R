library(data.table)
library(cowplot)
library(ggrepel)
# devtools::install_github('boxiangliu/manhattan')
library(manhattan)

smr_fn='../processed_data/finemap/smr/UKBB_HCASMC_smr_results.out.smr'
fm_fn='../processed_data/finemap/finemap_mike/UKBB_GWAS1KG_EXOME_CAD_SOFT_META_PublicRelease_300517_txt_gz_hcasmc_eqtls_txt_gz.txt'
fig_dir='../figures/finemap/locuszoom/manhattan/'
if(!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}

read_finemap=function(fm_fn){
	fm=fread(fm_fn)[,list(chrom=chrom,pos,gene_name=gene,y=clpp_score,method='eCAVIAR')]
	fm[,rank:=rank(-y,ties.method='random'),by='gene_name']
	fm=fm[rank==1]
	fm$rank=NULL
	return(fm)
}


smr=fread(smr_fn)[,list(chrom=ProbeChr,gene_name=Gene,pos=Probe_bp,y=log10(p_SMR),method='SMR')]
fm=read_finemap(fm_fn)


data=rbind(smr,fm)

data[,method:=factor(method,level=c('SMR','eCAVIAR'))]
data[,label:=ifelse( (method=='SMR'&y<log10(5e-5)) | (method=='eCAVIAR'&y>0.05),gene_name,'')]
data[,chrom:=paste0('chr',chrom)]

dummy=data.table(method=c('SMR','eCAVIAR'),y=c(log10(5e-5),0.05))

p=manhattan(data,build='hg19')+
	facet_grid(method~.,scale='free_y')+
	scale_y_continuous(labels=function(x){abs(x)})+
	geom_text(aes(label=label),hjust=-0.2)+
	geom_hline(data=dummy,aes(yintercept=y),color='red',linetype=2)+
	ylab(paste('-log10(P)                 CLPP'))
save_plot(sprintf('%s/manhattan.pdf',fig_dir),p,base_width=8,base_height=4)


