library(data.table)
library(cowplot)
library(ggrepel)
# devtools::install_github('boxiangliu/manhattan')
library(manhattan)
library(stringr)

smr_fn='../processed_data/finemap/smr/UKBB_HCASMC_eqtl_smr_results5e-05.out.smr'
fm_fn='../processed_data/finemap/finemap_mike/UKBB_GWAS1KG_EXOME_CAD_SOFT_META_PublicRelease_300517_txt_gz_finemap_clpp_status.txt'
fig_dir='../figures/finemap/locuszoom/manhattan/'
out_dir='../processed_data/finemap/locuszoom/manhattan/'
if(!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}
if(!dir.exists(out_dir)){dir.create(out_dir,recursive=TRUE)}

read_smr=function(smr_fn){
	fread(smr_fn)[,list(chrom=ProbeChr,gene_name=Gene,pos=Probe_bp,y=log10(p_SMR),gwas_logp=-log10(p_GWAS),eqtl_logp=-log10(p_eQTL),method='SMR')]
}

read_finemap=function(fm_fn,threshold=5e-5){
	if (grepl('clpp_status',fm_fn)) {
		tmp=fread(fm_fn,col.names=c('snp','eqtl_file','gwas_file','gene_name','conditional_level','n_tested_snps','clpp_score','gwas_logp','eqtl_logp'))
		tmp[,chrom:=str_split_fixed(snp,'_',2)[,1]]
		tmp[,pos:=as.integer(str_split_fixed(snp,'_',2)[,2])]
		fm=tmp[grepl('eqtl',eqtl_file),list(chrom,pos,gene_name,y=clpp_score,gwas_logp,eqtl_logp,method='eCAVIAR')]
	} else {
		fm=fread(fm_fn)[,list(chrom=chrom,pos,gene_name=gene,y=clpp_score,gwas_logp=gwas_log_pval,eqtl_logp=eqtl_log_pval,method='eCAVIAR')]
	}

	set.seed(42)
	fm[,rank:=rank(-y,ties.method='random'),by='gene_name']
	fm=fm[rank==1]
	fm$rank=NULL
	fm=fm[gwas_logp> -log10(threshold)|eqtl_logp> -log10(threshold)]
	return(fm)
}


smr=read_smr(smr_fn)
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
saveRDS(list(data,dummy,p),sprintf('%s/manhattan.rds',out_dir))
save_plot(sprintf('%s/manhattan.pdf',fig_dir),p,base_width=8,base_height=4)


