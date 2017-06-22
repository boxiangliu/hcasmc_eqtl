library(data.table) 
library(stringr)
library(cowplot)

rasqual_dir='../processed_data/finemap/finemap/rasqual_finemap/n_causal_max_1/'
gwas_dir='../processed_data/finemap/finemap/howson_finemap/n_causal_max_1/'

gwas_in_dir='../processed_data/finemap/finemap/howson_data/'
rasqual_in_dir='../processed_data/finemap/finemap/rasqual_data/input/'

fig_dir='../figures/finemap/finemap/howson_rasqual_coloc/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

container=list()
for (f in list.files(rasqual_dir,'snp')){
	print(f)

	rasqual_fn=sprintf('%s/%s',rasqual_dir,f)
	rasqual=fread(rasqual_fn)

	lociID=str_extract(f,'(?<=^)(chr.+?)(?=\\.)')
	gene_id=str_extract(f,'(?<=\\.)(.+?)(?=\\.snp)')

	gwas_fn=sprintf('%s/%s.snp',gwas_dir,lociID)
	gwas=fread(gwas_fn)

	merged=merge(gwas[,list(snp,snp_prob)],
		rasqual[,list(snp,snp_prob)],by='snp',
		suffixes=c('_gwas','_rasqual'))

	merged[,clpp:=snp_prob_gwas*snp_prob_rasqual]

	merged$gene_id=gene_id
	merged$lociID=lociID

	container[[f]]=merged
}

clpp=Reduce(rbind,container)

rclpp=clpp[,list(rclpp=sum(clpp)),by=c('gene_id','lociID')]
setorder(rclpp,-rclpp)

container=list()
for (i in 1:nrow(rclpp)){
	print(i)

	fid=rclpp[i,gene_id]
	id=rclpp[i,lociID]
	text=rclpp[i,rclpp]

	to_plot=clpp[gene_id==fid&lociID==id,]

	p=ggplot(to_plot,aes(x=snp_prob_gwas,y=snp_prob_rasqual))+
	geom_point()+scale_x_sqrt()+xlab('GWAS causal prob.')+
	ylab('RASQUAL causal prob.')+ggtitle(paste(id,fid,text,sep='\n'))

	container[[i]]=p
}

pdf(sprintf('%s/howson_rasqual_coloc.pdf',fig_dir),height=4,width=4)
for (i in seq_along(container)) {print(container[[i]])}
dev.off()


container=list()
for (i in 1:nrow(rclpp)){
	print(i)

	fid=rclpp[i,gene_id]
	id=rclpp[i,lociID]
	text=rclpp[i,rclpp]

	gwas=fread(sprintf('%s/%s.z',gwas_in_dir,id),col.names=c('snpID','zscore'))
	temp=data.table(gwas[,str_split_fixed(snpID,'_',5)[,1:4]])

	setnames(temp,c('chr','pos','ref','alt'))
	gwas=cbind(gwas,temp)
	gwas[,pos:=as.integer(pos)]

	gwas[,p:=2*pnorm(-abs(zscore))]
	gwas[,logp:=-log10(p)]

	rasqual=fread(sprintf('%s/%s.%s.z',rasqual_in_dir,id,fid),,col.names=c('snpID','zscore'))
	temp=data.table(rasqual[,str_split_fixed(snpID,'_',5)[,1:4]])

	setnames(temp,c('chr','pos','ref','alt'))
	rasqual=cbind(rasqual,temp)
	rasqual[,pos:=as.integer(pos)]

	rasqual[,p:=2*pnorm(-abs(zscore))]
	rasqual[,logp:=-log10(p)]

	chr=unique(gwas$chr)
	to_plot=rbind(gwas[,list(pos,logp,data='GWAS')],rasqual[,list(pos,logp,data='RASQUAL')])
	p=ggplot(to_plot,aes(pos,logp))+geom_point(size=2,alpha=0.5)+facet_grid('data~.',scales="free_y")+xlab(paste0('chr',chr))+ylab('-log10(P-value)')+ggtitle(paste(fid,id,sep='\n'))

	container[[i]]=p
}

pdf(sprintf('%s/howson_rasqual_locuszoom.pdf',fig_dir))
for (i in seq_along(container)) {print(container[[i]])}
dev.off()
