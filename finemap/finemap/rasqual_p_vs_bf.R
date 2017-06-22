library(data.table)
library(stringr)
library(cowplot)


eqtl_in_dir='../processed_data/finemap/finemap/rasqual_data/input/'
eqtl_out_dir='../processed_data/finemap/finemap/rasqual_finemap/n_causal_max_1'
fig_dir='../figures/finemap/finemap/rasqual_p_vs_bf/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

z_fn=list.files(eqtl_in_dir,'.z')

container=list()
for (f in z_fn){
	message('INFO - ', f)

	z=fread(sprintf('%s/%s',eqtl_in_dir,f),col.names=c('snp','zscore'))

	bf_fn=str_replace(f,'z','snp')
	bf=fread(sprintf('%s/%s',eqtl_out_dir,bf_fn))
	
	stopifnot(setequal(z$snp,bf$snp))

	merged=merge(z,bf,by='snp')
	merged[,p:=2*pnorm(-abs(zscore))]
	merged[,logp:=-log10(p)]
	
	config_fn=str_replace(f,'z','config')
	config=fread(sprintf('%s/%s',eqtl_out_dir,config_fn),nrow=1)
	causal=config[,unlist(str_split(config,','))]
	merged[,causal:=snp%in%causal]

	id=str_replace(f,'.z','')
	p1=ggplot(merged,aes(logp,snp_log10bf,color=causal,alpha=causal))+
		geom_point(size=2)+
		xlab('-log10(P-value)')+ylab('Log10(BF)')+
		scale_color_manual(values=c(`TRUE`='red',`FALSE`='black'),name='Causal')+
		scale_alpha_manual(values=c(`TRUE`=1,`FALSE`=0.05),guide='none')+
		ggtitle(id)

	p2=ggplot(merged,aes(logp,snp_prob,color=causal,alpha=causal))+
		geom_point(size=2)+
		xlab('-log10(P-value)')+ylab('Causal Prob.')+
		scale_color_manual(values=c(`TRUE`='red',`FALSE`='black'),name='Causal')+
		scale_alpha_manual(values=c(`TRUE`=1,`FALSE`=0.05),guide='none')+
		ggtitle(id)

	p=plot_grid(p1,p2)

	container[[id]]=p
}

pdf(sprintf('%s/p_vs_bf.n_causal_max_1.pdf',fig_dir),height=4,width=8)
for (i in seq_along(container)) {print(container[[i]])}
dev.off()



