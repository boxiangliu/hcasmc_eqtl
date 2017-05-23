library(data.table)
library(stringr)
library(cowplot)

# Variables:
in_dir='../processed_data/compare_fastqtl_and_rasqual/top_snp_per_gene/'
atac_fn='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks/HCASMC.merged.bed'
fig_dir='../figures/compare_fastqtl_and_rasqual/atacseq_overlap/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

# Read inputs: 
ras=fread(sprintf('%s/rasqual.txt',in_dir))
fas=fread(sprintf('%s/fastqtl.txt',in_dir))
atac=fread(atac_fn)


# Extend the peak boundaries by 500 each way: 
atac[,c('start','end'):=list(start-500,end+500)]


# Add start and end columns to RASQUAL:
ras[,c('start','end'):=pos]



# Add start and end columns to FastQTL:
tmp=as.data.table(fas[,str_split_fixed(sid,'_',5)[,1:4]])
setnames(tmp,c('chr','pos','ref','alt'))
tmp[,chr:=paste0('chr',chr)]
tmp[,pos:=as.numeric(pos)]
fas=cbind(fas,tmp)
fas[,c('start','end'):=pos]
rm(tmp)


# Overlap eQTL and ATACseq:
setkey(atac,chr,start,end)
tmp=unique(ras[,list(gene_id,pval)])
tmp[,gene_rank:=rank(pval)]
ras=merge(ras,tmp,by=c('gene_id','pval'))

tmp=unique(fas[,list(gene_id,pval)])
tmp[,gene_rank:=rank(pval)]
fas=merge(fas,tmp,by=c('gene_id','pval'))

y=data.frame()
for (i in seq(1e3,18e3,1e3)){
	fas_top=fas[gene_rank<=i,]
	ras_top=ras[gene_rank<=i,]

	setkey(ras_top,chr,start,end)
	ras_merge=foverlaps(atac,ras_top,nomatch=0)
	ras_snp=length(unique(ras_merge$sid))/length(unique(ras_top$sid))
	ras_gene=length(unique(ras_merge$gene_id))/length(unique(ras_top$gene_id))
	ras_pair=nrow(unique(ras_merge[,list(sid,gene_id)]))/nrow(unique(ras_top[,list(sid,gene_id)]))

	setkey(fas_top,chr,start,end)
	fas_merge=foverlaps(atac,fas_top,nomatch=0)
	fas_snp=length(unique(fas_merge$sid))/length(unique(fas_top$sid))
	fas_gene=length(unique(fas_merge$gene_id))/length(unique(fas_top$gene_id))
	fas_pair=nrow(unique(fas_merge[,list(sid,gene_id)]))/nrow(unique(fas_top[,list(sid,gene_id)]))

	y=rbind(y,data.frame(ras_gene,ras_snp,ras_pair,fas_gene,fas_snp,fas_pair,top=i))
}

to_plot=melt(y,id.vars='top',measure.vars=c('ras_gene','fas_gene'),variable.name='type',value.name='pct')
p1=ggplot(to_plot,aes(top,pct,color=type))+geom_point()+geom_line()+scale_color_discrete(name='Method',labels=c('RASQUAL','FastQTL'))+xlab('Top Genes')+ylab('Proportion overlap')+ggtitle('Gene')
to_plot=melt(y,id.vars='top',measure.vars=c('ras_snp','fas_snp'),variable.name='type',value.name='pct')
p2=ggplot(to_plot,aes(top,pct,color=type))+geom_point()+geom_line()+scale_color_discrete(name='Method',labels=c('RASQUAL','FastQTL'))+xlab('Top SNPs')+ylab('Proportion overlap')+ggtitle('SNP')
to_plot=melt(y,id.vars='top',measure.vars=c('ras_pair','fas_pair'),variable.name='type',value.name='pct')
p3=ggplot(to_plot,aes(top,pct,color=type))+geom_point()+geom_line()+scale_color_discrete(name='Method',labels=c('RASQUAL','FastQTL'))+xlab('Top Pairs')+ylab('Proportion overlap')+ggtitle('Association')


# Calculate the number of SNPs per gene:
y=data.frame()
for (i in seq(1e3,18e3,1e3)){
	x1=ras[gene_rank<=i&(i-1000)<=gene_rank,.N]
	x2=fas[gene_rank<=i&(i-1000)<=gene_rank,.N]
	tmp=rbind(data.frame(N=x1,Method='RASQUAL',top=i),data.frame(N=x2,Method='FastQTL',top=i))
	y=rbind(tmp,y)
}

p4=ggplot(y,aes(x=top,y=N,color=Method))+geom_point()+geom_line()+xlab('Top Genes')+ylab('Associated SNPs')
pdf(sprintf('%s/atacseq_overlap.pdf',fig_dir));p1;p2;p3;p4;dev.off()


# Calculate distance between eQTL and nearest ATACseq peak: 
y=atac[chr=='chr1',]
y[,summit:=(start+end)/2]

x=ras[chr=='chr1',]
o=outer(x$pos,y$summit,`-`)
m=apply(o,1,function(x) {x[abs(x)==min(abs(x))]})
m=Reduce(c,m)


x2=fas[chr=='chr1',]
o=outer(x2$pos,y$summit,`-`)
m2=apply(o,1,function(x) {x[abs(x)==min(abs(x))]})
m2=Reduce(c,m2)

bins=cut(m,breaks=seq(-1e4,1e4,2e4/40))
table1=table(bins)
table1=table1/sum(table1)
bins=cut(m2,breaks=seq(-1e4,1e4,2e4/40))
table2=table(bins)
table2=table2/sum(table2)

to_plot=rbind(data.frame(as.data.frame(table1),Method='RASQUAL'),data.frame(as.data.frame(table2),Method='FastQTL'))
p5=ggplot(to_plot,aes(bins,Freq,fill=Method))+geom_bar(stat='identity',position=position_dodge(width=0.8),alpha=0.95)+scale_x_discrete(name='Distance',labels=seq(-1e4,1e4,2e4/40))+theme(axis.text.x=element_text(angle=90))
pdf(sprintf('%s/dist_to_peaks.pdf',fig_dir));p5;dev.off()





















