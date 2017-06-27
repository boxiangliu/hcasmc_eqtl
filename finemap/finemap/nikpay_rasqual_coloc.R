library(data.table)
library(stringr)
library(cowplot)

gwas_fn='../data/gwas/CARDIoGRAMplusC4D/cad.add.160614.website.txt'
gencode_fn='../data/gtex/gencode.v19.genes.v6p.hg19.bed'
eqtl_dir='../processed_data/rasqual/output_pval/'
ld_block_fn='/users/bliu2/tools/LDetect/ldetect-data/EUR/fourier_ls-all.bed'


fig_dir='../figures/finemap/finemap/nikpay_rasqual_coloc/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}



gwas=fread(gwas_fn,select=c(1:6,8:11),
	col.names=c('rsid','chr','pos','effect_allele',
		'other_allele','freq','model','beta','se','p'))
gwas[,c('start','end'):=as.integer(pos)]
gwas[,chr:=paste0('chr',chr)]

ld_block=fread(ld_block_fn,header=TRUE,
	col.names=c('chr','start','end'))
ld_block[,start:=start+1]
ld_block[,lociID:=sprintf('%s_%s_%s',chr,start,end)]

setkey(gwas,chr,start,end)
setkey(ld_block,chr,start,end)
overlap=foverlaps(gwas[,list(chr,start,end,pos,rsid,p)],
	ld_block,nomatch=0)
stopifnot(nrow(overlap)==nrow(gwas))

overlap[,c('i.start','i.end'):=NULL]
overlap[,logp:=-log10(p)]

min_p=overlap[,list(min_p=min(p)),by='lociID']
sig_loci=min_p[min_p<1e-5,lociID]
sig_loci=data.table(sig_loci,str_split_fixed(sig_loci,'_',3))
setnames(sig_loci,c('lociID','chr','start','end'))
sig_loci[,c('start','end'):=list(as.integer(start),as.integer(end))]

gencode=fread(gencode_fn,col.names=c('chr','start','end',
	'strand','gene_id','gene_name','type'))

window_size=1e6
gencode[,window_start:=ifelse(strand=='+',start-window_size,end-window_size)]
gencode[,window_end:=ifelse(strand=='+',start+window_size,end+window_size)]

setkey(gencode,chr,window_start,window_end)
setkey(sig_loci,chr,start,end)
overlap1=foverlaps(sig_loci,gencode,nomatch=0)
overlap1=overlap1[type%in%c('protein_coding','lincRNA')]

eqtl_threshold=1e-6
container=list()
n=0
for (i in 1:nrow(overlap1)){
	message(i)
	chr=overlap1[i,chr]
	gene_id=overlap1[i,gene_id]
	gene_name=overlap1[i,gene_name]

	eqtl_fn=sprintf('%s/%s/%s_%s.pval.txt',eqtl_dir,chr,gene_id,gene_name)
	if (file.exists(eqtl_fn)){
		eqtl=fread(eqtl_fn)
	} else {
		warning(eqtl_fn,' does not exists!')
	}

	if (min(eqtl$pval)<eqtl_threshold){
		n=n+1
		container[[n]]=eqtl
	}
}

eqtl=Reduce(rbind,container)
eqtl=unique(eqtl)
eqtl[,c('start','end'):=as.integer(pos)]
setkey(eqtl,chr,start,end)
eqtl_overlap=foverlaps(eqtl,sig_loci,nomatch=0)
eqtl_overlap[,c('i.start','i.end','start','end'):=NULL]
eqtl_overlap[,logp:=-log10(pval)]

iterator=unique(eqtl_overlap[,list(lociID,fid)])

setkey(eqtl_overlap,lociID,fid)
container=list()
for (i in 1:nrow(iterator)){
	id=iterator[i,lociID]
	gene_id=iterator[i,fid]
	rasqual=eqtl_overlap[lociID==id&fid==gene_id,list(chr,pos,logp,data='RASQUAL')]

	gwas=overlap[lociID==id,list(chr,pos,logp,data='GWAS')]

	to_plot=rbind(rasqual,gwas)

	p=ggplot(to_plot,aes(pos,logp))+geom_point(size=2,alpha=0.2)+
		facet_grid('data~.',scales="free_y")+
		xlab(paste0(chr))+ylab('-log10(P-value)')+
		ggtitle(paste(gene_id,id,sep='\n'))

	container[[i]]=p
}

pdf(sprintf('%s/nikpay_rasqual_locuszoom.pdf',fig_dir))
for (i in seq_along(container)) {print(container[[i]])}
dev.off()
