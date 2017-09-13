library(data.table)
library(stringr)
library(cowplot)
library(rtracklayer)
library(ggrepel)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)

atac_dir='../data/atacseq/fbs/'

bigwig_fn=c("../data/atacseq/fbs/1346/out/signal/macs2/rep1/CA1346_L2_TCCTGAGC_L002_R1.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig",
	"../data/atacseq/fbs/1508/out/signal/macs2/pooled_rep/CA1508_L1_TAAGGCGA_L001_R1.trim.PE2SE.nodup.tn5_pooled.pf.fc.signal.bigwig",
	"../data/atacseq/fbs/1522/out/signal/macs2/rep1/CA1522_S8_concat_R1_001.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig",
	"../data/atacseq/fbs/200212/out/signal/macs2/rep1/200212_L2_CGTACTAG_L002_R1.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig",
	"../data/atacseq/fbs/2108/out/signal/macs2/rep1/CA2108_S4_concat_R1_001.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig",
	"../data/atacseq/fbs/2305/out/signal/macs2/pooled_rep/CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.nodup.tn5_pooled.pf.fc.signal.bigwig",
	"../data/atacseq/fbs/2356/out/signal/macs2/rep1/CA2356_S7_concat_R1_001.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig",
	"../data/atacseq/fbs/2510/out/signal/macs2/rep1/CA2510_S5_concat_R1_001.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig",
	"../data/atacseq/fbs/2989/out/signal/macs2/rep1/CA2989_S7_concat_R1_001.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig")

chromHMM_fn='../data/chromatin/chromHMM/out/HCASMC_10_dense.bed'
dbsnp_fn='../data/dbsnp/snp147.txt'

fig_dir='../figures/finemap/finemap/atacseq_overlap/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

# Function: 
get_bw=function(f,seq,start,end){
	gr=GRanges(seqnames = seq,
		ranges = IRanges(start = start, end = end))

	bw=import(f,which=gr)
	bw=as.data.frame(bw)
	setDT(bw)
	
	temp1=bw[,list(seqnames,start,score)]
	temp2=bw[,list(seqnames,end,score)]
	setnames(temp1,'start','pos')
	setnames(temp2,'end','pos')
	bw2=rbind(temp1,temp2)

	setorder(bw2,seqnames,pos)
	return(bw2)
}


chromHMM2GRanges=function(file){
	chromHMM=fread(file,sep='\t',skip=1,
		col.names=c('chr','start','end','name','score','strand','thickStart','thickEnd','itemRgb'))
	chromHMM=GenomicRanges::makeGRangesFromDataFrame(chromHMM,keep.extra.columns=TRUE)
	return(chromHMM)
}


snp2GRanges=function(file){
	snp=fread(file,sep='\t',skip=1,select=c(2,3,4,5),
		col.names=c('chr','start','end','name'))
	snp=GenomicRanges::makeGRangesFromDataFrame(snp,keep.extra.columns=TRUE)
	return(snp)
}


combine_data=function(bigwig_fn,nikpay,howson,chromHMM_fn,dbsnp_fn){
	gen='hg19'

	tracks=list()
	for (f in bigwig_fn){
		print(f)
		sample=str_extract(f,'(?<=/)([0-9]+?)(?=/)')

		bw_track = DataTrack(range = f, genome = gen, 
			window = -1, name = paste0('ATAC ',sample), type='histogram',frame=TRUE,col.frame='black')
	
		tracks[[sample]]=bw_track
	}

	print('howson')
	# ho=howson
	# ho_track = DataTrack(data = ho$logp, chromosome = ho$chr, start = ho$pos, 
	# 	end = ho$pos, id = ho$rsid, genome = gen, name = 'Howson', type='p',frame=TRUE,col.frame='black')
	tracks[['Howson']] = ho_track


	print('nikpay')
	# ni=nikpay
	# ni_track = DataTrack(data = ni$logp, chromosome = ni$chr, start = ni$pos, 
	# 	end = ni$pos, id = ni$rsid, genome = gen, name = 'Nikpay', type='p',frame=TRUE,col.frame='black')
	tracks[['Nikpay']] = ni_track


	print('chromHMM')
	chromHMM=chromHMM2GRanges(chromHMM_fn)
	ch_track = AnnotationTrack(range = chromHMM, genome = gen, 
		name = 'ChromHMM', feature = chromHMM$itemRgb, id = chromHMM$name, 
		stacking = 'dense', collapse = FALSE, frame=TRUE,col.frame='black',groupAnnotation = 'id')

	feat = unique(feature(ch_track))
	featCol = setNames(as.list(rgb(t(sapply(strsplit(feat, ","),
		as.numeric)), maxColorValue=255)), feat)
	displayPars(ch_track) = featCol

	tracks[['chromHMM']] = ch_track

	print('Gene model')
	# grtrack = GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene, genome = gen,name = 'Gene Model',  transcriptAnnotation = "symbol",frame=TRUE,col.frame='black')
	# symbols = unlist(mapIds(org.Hs.eg.db, gene(grtrack), "SYMBOL", "ENTREZID", multiVals = "first"))
	# symbol(grtrack) = symbols[gene(grtrack)]
	tracks[['gene_model']] = grtrack


	print('dbSNP')
	snp=snp2GRanges(dbsnp_fn)
	snp_track = AnnotationTrack(range = snp, genome = gen, name = 'SNPs', id = snp$name, 
		stacking = 'full', collapse = FALSE, frame=TRUE,col.frame='black', groupAnnotation = 'id')
	tracks[['snp']] = snp_track


	print('Genomic location')
	gtrack = GenomeAxisTrack(labelPos = "below")
	tracks[['genomic_location']] = gtrack

	return(tracks)
}

# grtrack=tracks[['gene_model']]
# ho_track=tracks[['Howson']]
# ni_track=tracks[['Nikpay']]

plot_ase=function(x,title){
	pi=x[,pi]
	freq=x[,freq]
	ref=x[,ref]
	alt=x[,alt]

	to_plot=data.frame(allele=c(sprintf('ref: %s (%.02f)',ref,1-freq),sprintf('alt: %s (%.02f)',alt,freq)),ase=c(1-pi,pi))
	p=ggplot(to_plot,aes(allele,ase,label=sprintf('%.02f%%',ase*100)))+geom_bar(stat='identity')+geom_text(nudge_y=0.02)+ggtitle(title)
	return(p)
}


# Read Nikpay data:
nikpay_fn='../data/gwas/CARDIoGRAMplusC4D/cad.add.160614.website.txt'
nikpay=fread(nikpay_fn,select=c(1:6,8:11),
	col.names=c('rsid','chr','pos','effect_allele',
		'other_allele','freq','model','beta','se','p'))
nikpay[,chr:=paste0('chr',chr)]
nikpay[,logp:=-log10(p)]


# Read Howson data:
howson_fn='../data/gwas/howson_2017/Howson-JMM_CHD_Mixed_2017.norm.in1kgp3.txt'
howson=fread(howson_fn,header=TRUE)
temp=howson[,str_split_fixed(chrpos_b37,':',2)]
howson[,c('chr','pos'):=list(temp[,1],as.integer(temp[,2]))]
howson[,logp:=-log10(p)]



# Combine all data:
tracks=combine_data(bigwig_fn,nikpay,howson,chromHMM_fn,dbsnp_fn)



# TCF21:
chr='chr6'
start=134e6
end=134.4e6
eqtl_fn='../processed_data/rasqual/output_pval/chr6/ENSG00000118526.6_TCF21.pval.txt'


eqtl=fread(eqtl_fn)
eqtl[,logp:=-log10(pval)]

eqtl_track = DataTrack(data = eqtl$logp, start = eqtl$pos, 
	end = eqtl$pos, chromosome = chr, genome = 'hg19', name = 'RAS', type='p',frame=TRUE,col.frame='black')

tracks[['RASQUAL']]=eqtl_track
tracks=tracks[c(1:11,16,12:15)]

pdf(sprintf('%s/tcf21.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(tracks[-15], chromosome = chr, from = start, to = end)
dev.off()


start=134200800
end=134215600
ht = HighlightTrack(trackList=tracks, start = c(134209812,134214500), width = 50, chromosome= chr)
pdf(sprintf('%s/tcf21.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(list(ht), chromosome = chr, from = start, to = end,showFeatureId=FALSE,showId=TRUE)
dev.off()


# FES:
chr='chr15'
start=905e5
end=920e5
eqtl_fn='../processed_data/rasqual/output_pval/chr15/ENSG00000182511.7_FES.pval.txt'

eqtl=fread(eqtl_fn)
eqtl[,logp:=-log10(pval)]

eqtl_track = DataTrack(data = eqtl$logp, start = eqtl$pos, 
	end = eqtl$pos, chromosome = chr, genome = 'hg19', name = 'RAS', type='p',frame=TRUE,col.frame='black')

tracks[['RASQUAL']]=eqtl_track

pdf(sprintf('%s/fes.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(tracks[-15], chromosome = chr, from = start, to = end)
dev.off()


start=913e5
end=915e5
pdf(sprintf('%s/fes.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(tracks[-15], chromosome = chr, from = start, to = end)
dev.off()


start=91400300
end=91448950
pdf(sprintf('%s/fes.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(tracks[-15], chromosome = chr, from = start, to = end)
dev.off()

start=91.427e6
end=91.43e6
pos = c(91427612,91427872,91428197,91428290,
	91428522,91428589,91428636,91428955,
	91429042,91429176,91429198,91429287)
ht = HighlightTrack(trackList=tracks, start = pos-5, width = 10, chromosome= chr)
pdf(sprintf('%s/fes.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(ht, chromosome = chr, from = start, to = end)
dev.off()

eqtl[rsid%in%c('rs35346340','rs2071382','rs11539637','rs11330240','rs7183988','rs7177338','rs1894400','rs1894401','rs7497304','rs34029266','rs4932373'),list(pval)]
nikpay[rsid%in%c('rs12906125','rs35346340','rs2071382','rs11539637','rs11330240','rs7183988','rs7177338','rs1894400','rs1894401','rs7497304','rs34029266','rs4932373'),list(p)]
nikpay[chr=='chr15'&pos==91429196,list(rsid,effect_allele,other_allele,p)]


# SIPA1:
chr='chr11'
start=65e6
end=66e6
eqtl_fn='../processed_data/rasqual/output_pval/chr11/ENSG00000213445.4_SIPA1.pval.txt'
eqtl=fread(eqtl_fn)
eqtl[,logp:=-log10(pval)]

eqtl_track = DataTrack(data = eqtl$logp, start = eqtl$pos, 
	end = eqtl$pos, chromosome = chr, genome = 'hg19', name = 'RAS', type='p',frame=TRUE,col.frame='black')

tracks[['RASQUAL']]=eqtl_track

pdf(sprintf('%s/sipa1.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(tracks[-15], chromosome = chr, from = start, to = end)
dev.off()


start=65375e3
end=65475e3

pdf(sprintf('%s/sipa1.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(tracks[-15], chromosome = chr, from = start, to = end)
dev.off()


start=65400000
end=65412500
pos = c(65405600,65406869,65408937)
ht = HighlightTrack(trackList=tracks, start = pos-5, width = 10, chromosome= chr)
pdf(sprintf('%s/sipa1.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(ht, chromosome = chr, from = start, to = end)
dev.off()



start=65380000
end=65400000
pos=65391317
ht = HighlightTrack(trackList=tracks, start = pos-5, width = 10, chromosome= chr)
pdf(sprintf('%s/sipa1.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(ht, chromosome = chr, from = start, to = end)
dev.off()



# MRAS:
chr='chr3'
start=1.375e8
end=1.385e8
eqtl_fn='../processed_data/rasqual/output_pval/chr3/ENSG00000158186.8_MRAS.pval.txt'
eqtl=fread(eqtl_fn)
eqtl[,logp:=-log10(pval)]


eqtl_track = DataTrack(data = eqtl$logp, start = eqtl$pos, 
	end = eqtl$pos, chromosome = chr, genome = 'hg19', name = 'RAS', type='p',frame=TRUE,col.frame='black')

tracks[['RASQUAL']]=eqtl_track

pdf(sprintf('%s/mras.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(tracks[-15], chromosome = chr, from = start, to = end)
dev.off()


start=138000000
end=138250000
pdf(sprintf('%s/mras.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(tracks[-15], chromosome = chr, from = start, to = end)
dev.off()


start=138.05e6
end=138.125e6
pdf(sprintf('%s/mras.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(tracks, chromosome = chr, from = start, to = end)
dev.off()

start=138.05e6
end=138.06e6
pos=138052754
ht = HighlightTrack(trackList=tracks, start = pos-5, width = 10, chromosome= chr)
pdf(sprintf('%s/mras.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(ht, chromosome = chr, from = start, to = end)
dev.off()

start=138.07e6
end=138.075e6
pos=138070901
ht = HighlightTrack(trackList=tracks, start = pos-5, width = 10, chromosome= chr)
pdf(sprintf('%s/mras.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(ht, chromosome = chr, from = start, to = end)
dev.off()


start=138.085e6
end=138.102e6
pos=c(138087467, 138088064, 138092889, 138095525, 138096097, 138099161)
ht = HighlightTrack(trackList=tracks, start = pos-5, width = 10, chromosome= chr)
pdf(sprintf('%s/mras.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(ht, chromosome = chr, from = start, to = end)
dev.off()


start=138.107e6
end=138.123e6
pos=c(138108352, 138111751, 138119952, 138121920, 138122122)
ht = HighlightTrack(trackList=tracks, start = pos-5, width = 10, chromosome= chr)
pdf(sprintf('%s/mras.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(ht, chromosome = chr, from = start, to = end)
dev.off()


# SNHG18:
chr='chr5'
eqtl_fn='../processed_data/rasqual/output_pval/chr5/ENSG00000250786.1_SNHG18.pval.txt'
eqtl=fread(eqtl_fn)
eqtl[,logp:=-log10(pval)]


eqtl_track = DataTrack(data = eqtl$logp, start = eqtl$pos, 
	end = eqtl$pos, chromosome = chr, genome = 'hg19', name = 'RAS', type='p',frame=TRUE,col.frame='black')

tracks[['RASQUAL']]=eqtl_track

start=95e5
end=96e5
pdf(sprintf('%s/snhg18.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(tracks[-15], chromosome = chr, from = start, to = end)
dev.off()


start=9.55e6
end=9.56e6
pos=9556694
ht = HighlightTrack(trackList=tracks, start = pos-5, width = 10, chromosome= chr)
pdf(sprintf('%s/snhg18.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(ht, chromosome = chr, from = start, to = end)
dev.off()


start=9.543e6
end=9.545e6
pos=9543580
ht = HighlightTrack(trackList=tracks, start = pos-5, width = 10, chromosome= chr)
pdf(sprintf('%s/snhg18.%s_%s_%s.pdf',fig_dir,chr,start,end),height=10,width=8)
plotTracks(ht, chromosome = chr, from = start, to = end)
dev.off()

eqtl[rsid=='rs1508798',list(paste(chr,pos,ref,alt,sep=':'))]




