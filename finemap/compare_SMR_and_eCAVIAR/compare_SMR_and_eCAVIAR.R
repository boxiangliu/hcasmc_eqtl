library(data.table)
library(cowplot)
library(ggrepel)
# devtools::install_github('boxiangliu/manhattan')
library(manhattan)
library(stringr)

smr_fn='../processed_data/finemap/smr/UKBB_HCASMC_eqtl_smr_results5e-05.out.smr'
fm_fn='../processed_data/finemap/finemap_mike/UKBB_GWAS1KG_EXOME_CAD_SOFT_META_PublicRelease_300517_txt_gz_finemap_clpp_status.txt'
fig_dir='../figures/finemap/compare_SMR_and_eCAVIAR/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir)}

if(!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}
if(!dir.exists(out_dir)){dir.create(out_dir,recursive=TRUE)}

read_smr=function(smr_fn){
	smr = fread(smr_fn)[,list(chrom=ProbeChr,gene_name=Gene,pos=Probe_bp,y=-log10(p_SMR),gwas_logp=-log10(p_GWAS),eqtl_logp=-log10(p_eQTL),method='SMR')]
	smr = smr[!is.na(y)]
	return(smr)
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

merged = merge(smr[,list(gene_name,smr = y)],fm[,list(gene_name, fm = y)],by='gene_name')
merged[,smr_pval := 10^-smr]
merged[,smr_bin := cut(smr_pval,breaks = 10^-c(0,1,2,Inf))]
merged[,smr_bin := factor(smr_bin, rev(levels(smr_bin)))]
p = ggplot(merged,aes(smr_bin,fm)) + 
	geom_boxplot() + 
	scale_y_log10() + 
	scale_x_discrete(breaks = c("(0.1,1]","(0.01,0.1]","(0,0.01]"), 
		labels = c(expression('0.1 < P'<='1'), 
			expression('0.01 < P'<='0.1'), 
			expression('0 < P'<='0.01'))) + 
	xlab('SMR P-values') + 
	ylab('eCAVIAR colocalization probability')
fig_fn = sprintf('%s/compare_SMR_and_eCAVIAR.pdf',fig_dir)
save_plot(fig_fn, p)