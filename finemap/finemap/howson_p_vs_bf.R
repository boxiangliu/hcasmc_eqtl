library(data.table)
library(stringr)
library(cowplot)
library(foreach)
library(doMC)
registerDoMC(cores=20)


gwas_fn='../data/gwas/howson_2017/Howson-JMM_CHD_Mixed_2017.norm.in1kgp3.txt'
ld_block_fn='/users/bliu2/tools/LDetect/ldetect-data/EUR/fourier_ls-all.bed'
fig_dir='../figures/finemap/finemap/p_vs_bf/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

gwas=fread(gwas_fn,header=TRUE)

temp=gwas[,str_split_fixed(chrpos_b37,':',2)]
gwas[,c('chr','pos'):=list(temp[,1],temp[,2])]
gwas[,c('start','end'):=as.integer(pos)]

gwas[,zscore:=beta/se]

ld_block=fread(ld_block_fn,header=TRUE,col.names=c('chr','start','end'))
ld_block[,start:=start+1]
ld_block[,lociID:=sprintf('%s_%s_%s',chr,start,end)]

setkey(gwas,chr,start,end)
setkey(ld_block,chr,start,end)
overlap=foverlaps(gwas[,list(chr,start,end,rsid,snpID,
	ref,alt,p,n,chrpos_b37,zscore)],
	ld_block,nomatch=0)

overlap[,pos:=i.start]
overlap[,c('i.start','i.end'):=NULL]

min_p=overlap[,list(min_p=min(p)),by='lociID']
sig_loci=min_p[min_p<1e-6,lociID]


in_dir='../processed_data/finemap/finemap/finemap/default/'
container=list()
for (id in sig_loci){
	message('INFO - ', id)
	bf=fread(sprintf('%s/%s.snp',in_dir,id))

	merged=merge(overlap[lociID==id,list(lociID,chrpos_b37,p)],
		bf,by.x='chrpos_b37',by.y='snp')
	merged[,logp:=-log10(p)]

	config=fread(sprintf('%s/%s.config',in_dir,id),nrow=1)
	causal=config[,unlist(str_split(config,','))]
	
	merged[,causal:=chrpos_b37%in%causal]
	
	p1=ggplot(merged,aes(logp,snp_log10bf,color=causal))+
		geom_point(size=2,alpha=0.5)+
		xlab('-log10(P-value)')+ylab('Log10(BF)')+
		scale_color_manual(values=c(`TRUE`='red',`FALSE`='black'),
		name='Causal')+ggtitle(id)

	p2=ggplot(merged,aes(logp,snp_prob,color=causal))+
		geom_point(size=2,alpha=0.5)+
		xlab('-log10(P-value)')+ylab('Causal Prob.')+
		scale_color_manual(values=c(`TRUE`='red',`FALSE`='black'),
		name='Causal')+ggtitle(id)

	p=plot_grid(p1,p2)

	container[[id]]=p
}

pdf(sprintf('%s/p_vs_bf.pdf',fig_dir),height=4,width=8)
for (i in seq_along(container)) {print(container[[i]])}
dev.off()


in_dir='../processed_data/finemap/finemap/finemap//n_causal_max_1/'
container=list()
for (id in sig_loci){
	message('INFO - ', id)
	bf=fread(sprintf('%s/%s.snp',in_dir,id))
	
	merged=merge(overlap[lociID==id,list(lociID,chrpos_b37,p)],
		bf,by.x='chrpos_b37',by.y='snp')
	merged[,logp:=-log10(p)]
	
	config=fread(sprintf('%s/%s.config',in_dir,id),nrow=1)
	causal=config[,unlist(str_split(config,','))]
	
	merged[,causal:=chrpos_b37%in%causal]
	
	p1=ggplot(merged,aes(logp,snp_log10bf,color=causal))+
		geom_point(size=2,alpha=0.5)+
		xlab('-log10(P-value)')+ylab('Log10(BF)')+
		scale_color_manual(values=c(`TRUE`='red',`FALSE`='black'),
		name='Causal')+ggtitle(id)

	p2=ggplot(merged,aes(logp,snp_prob,color=causal))+
		geom_point(size=2,alpha=0.5)+
		xlab('-log10(P-value)')+ylab('Causal Prob.')+
		scale_color_manual(values=c(`TRUE`='red',`FALSE`='black'),
		name='Causal')+ggtitle(id)

	p=plot_grid(p1,p2)

	container[[id]]=p
}

pdf(sprintf('%s/p_vs_bf.n_causal_max_1.pdf',fig_dir),height=4,width=8)
for (i in seq_along(container)) {print(container[[i]])}
dev.off()



in_dir='../processed_data/finemap/finemap/finemap//n_causal_max_2/'
container=list()
for (id in sig_loci){
	message('INFO - ', id)
	bf=fread(sprintf('%s/%s.snp',in_dir,id))
	
	merged=merge(overlap[lociID==id,list(lociID,chrpos_b37,p)],
		bf,by.x='chrpos_b37',by.y='snp')
	merged[,logp:=-log10(p)]
	
	config=fread(sprintf('%s/%s.config',in_dir,id),nrow=1)
	causal=config[,unlist(str_split(config,','))]
	
	merged[,causal:=chrpos_b37%in%causal]
	
	p1=ggplot(merged,aes(logp,snp_log10bf,color=causal))+
		geom_point(size=2,alpha=0.5)+
		xlab('-log10(P-value)')+ylab('Log10(BF)')+
		scale_color_manual(values=c(`TRUE`='red',`FALSE`='black'),
		name='Causal')+ggtitle(id)

	p2=ggplot(merged,aes(logp,snp_prob,color=causal))+
		geom_point(size=2,alpha=0.5)+
		xlab('-log10(P-value)')+ylab('Causal Prob.')+
		scale_color_manual(values=c(`TRUE`='red',`FALSE`='black'),
		name='Causal')+ggtitle(id)

	p=plot_grid(p1,p2)

	container[[id]]=p
}

pdf(sprintf('%s/p_vs_bf.n_causal_max_2.pdf',fig_dir),height=4,width=8)
for (i in seq_along(container)) {print(container[[i]])}
dev.off()
