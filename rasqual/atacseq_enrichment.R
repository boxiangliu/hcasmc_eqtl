library(data.table)
library(gap)

eqtl_fn='../processed_data/rasqual/output_merged/expressed_genes.pval.txt'
atacseq_fn='../data/atacseq/fbs/2305/out/peak/macs2/overlap/CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.nodup.tn5_pooled.pf.pval0.1.500K.naive_overlap.narrowPeak.gz'
fig_dir='../figures/rasqual/atacseq_enrichment/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir)}

read_eqtl=function(fn){
	eqtl=fread(fn,
		select=c(1,2,3,4,5,6,26,27,28),
		col.names=c('gene_name','rsid','chr','pos','ref','alt','pval','rank','gene_id'))
	return(eqtl)
}

get_top_eqtl=function(eqtl){
	top_eqtl=eqtl[rank==1]
	return(top_eqtl)
}

read_atacseq=function(fn){
	atacseq=fread(sprintf('zcat %s',fn),select=c(1,2,3),col.names=c('chr','start','end'))
	return(unique(atacseq))
}

overlap_eqtl_and_atacseq=function(top_eqtl,atacseq){
	top_eqtl[,c('start','end'):=list(pos,pos)]
	setkey(atacseq,chr,start,end)
	overlap=foverlaps(top_eqtl,atacseq,nomatch=0)
	return(overlap)
}

make_qqplot=function(top_eqtl,overlap,threshold = 1e-16){
	res1=qqunif(top_eqtl$pval)
	res2=qqunif(overlap$pval)
	plot(res1,col='blue',pch=16,cex=1.5,xlab='-log10(Expected)',ylab='-log10(Observed)')
	points(res2,col='red',pch=16,cex=1.5)
	abline(a=0,b=1,col='red')
	legend('topleft',legend=c('eQTL','eQTL in ATAC-seq'),pch=16,col=c('blue','red'))
}

eqtl=read_eqtl(eqtl_fn)
top_eqtl=get_top_eqtl(eqtl)
atacseq=read_atacseq(atacseq_fn)
overlap=overlap_eqtl_and_atacseq(top_eqtl,atacseq)
pval = wilcox.test(top_eqtl$pval,overlap$pval)$p.value # 9.221e-05
label = sprintf('Wilcoxon Rank-Sum Test P-value = %.02e',pval)

pdf(sprintf('%s/atacseq_enrichment.pdf',fig_dir))
make_qqplot(top_eqtl,overlap)
text(x=2,y=150,label=label)
dev.off()