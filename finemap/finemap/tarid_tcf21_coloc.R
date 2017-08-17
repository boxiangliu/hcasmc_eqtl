library(data.table)
library(Gviz)
library(GenomicRanges)
library(stringr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

# Variables:
bigwig_fn=c("../data/atacseq/fbs/1346/out/signal/macs2/rep1/CA1346_L2_TCCTGAGC_L002_R1.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig",
	"../data/atacseq/fbs/1508/out/signal/macs2/pooled_rep/CA1508_L1_TAAGGCGA_L001_R1.trim.PE2SE.nodup.tn5_pooled.pf.fc.signal.bigwig",
	"../data/atacseq/fbs/1522/out/signal/macs2/rep1/CA1522_S8_concat_R1_001.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig",
	"../data/atacseq/fbs/200212/out/signal/macs2/rep1/200212_L2_CGTACTAG_L002_R1.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig",
	"../data/atacseq/fbs/2108/out/signal/macs2/rep1/CA2108_S4_concat_R1_001.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig",
	"../data/atacseq/fbs/2305/out/signal/macs2/pooled_rep/CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.nodup.tn5_pooled.pf.fc.signal.bigwig",
	"../data/atacseq/fbs/2356/out/signal/macs2/rep1/CA2356_S7_concat_R1_001.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig",
	"../data/atacseq/fbs/2510/out/signal/macs2/rep1/CA2510_S5_concat_R1_001.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig",
	"../data/atacseq/fbs/2989/out/signal/macs2/rep1/CA2989_S7_concat_R1_001.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig")

fig_dir='../figures/finemap/finemap/tarid_tcf21_coloc/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

# Functions: 
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
	ho=howson
	ho_track = DataTrack(data = ho$logp, chromosome = ho$chr, start = ho$pos, 
		end = ho$pos, id = ho$rsid, genome = gen, name = 'Howson', type='p',frame=TRUE,col.frame='black')
	tracks[['Howson']] = ho_track


	print('nikpay')
	ni=nikpay
	ni_track = DataTrack(data = ni$logp, chromosome = ni$chr, start = ni$pos, 
		end = ni$pos, id = ni$rsid, genome = gen, name = 'Nikpay', type='p',frame=TRUE,col.frame='black')
	tracks[['Nikpay']] = ni_track


	# print('chromHMM')
	# chromHMM=chromHMM2GRanges(chromHMM_fn)
	# ch_track = AnnotationTrack(range = chromHMM, genome = gen, 
	# 	name = 'ChromHMM', feature = chromHMM$itemRgb, id = chromHMM$name, 
	# 	stacking = 'dense', collapse = FALSE, frame=TRUE,col.frame='black',groupAnnotation = 'id')

	# feat = unique(feature(ch_track))
	# featCol = setNames(as.list(rgb(t(sapply(strsplit(feat, ","),
	# 	as.numeric)), maxColorValue=255)), feat)
	# displayPars(ch_track) = featCol

	# tracks[['chromHMM']] = ch_track

	print('Gene model')
	grtrack = GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene, genome = gen,name = 'Gene Model',  transcriptAnnotation = "symbol",frame=TRUE,col.frame='black')
	symbols = unlist(mapIds(org.Hs.eg.db, gene(grtrack), "SYMBOL", "ENTREZID", multiVals = "first"))
	symbol(grtrack) = symbols[gene(grtrack)]
	tracks[['gene_model']] = grtrack


	# print('dbSNP')
	# snp=snp2GRanges(dbsnp_fn)
	# snp_track = AnnotationTrack(range = snp, genome = gen, name = 'SNPs', id = snp$name, 
	# 	stacking = 'full', collapse = FALSE, frame=TRUE,col.frame='black', groupAnnotation = 'id')
	# tracks[['snp']] = snp_track


	print('Genomic location')
	gtrack = GenomeAxisTrack(labelPos = "below")
	tracks[['genomic_location']] = gtrack

	return(tracks)
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

tracks=combine_data(bigwig_fn,nikpay,howson,-1,-1)


# TCF21:
chr='chr6'
start=1332e5
end=1352e5
eqtl_fn='../processed_data/rasqual/output_pval/chr6/ENSG00000118526.6_TCF21.pval.txt'

eqtl=fread(eqtl_fn)
eqtl[,logp:=-log10(pval)]


eqtl_track = DataTrack(data = eqtl$logp, start = eqtl$pos, 
	end = eqtl$pos, chromosome = chr, genome = 'hg19', name = 'TCF21', type='p',frame=TRUE,col.frame='black')

tracks[['TCF21']]=eqtl_track

# TARID:
eqtl_fn='../processed_data/rasqual/output_pval/chr6/ENSG00000227954.2_RP3-323P13.2.pval.txt'

eqtl=fread(eqtl_fn)
eqtl[,logp:=-log10(pval)]

eqtl_track = DataTrack(data = eqtl$logp, start = eqtl$pos, 
	end = eqtl$pos, chromosome = chr, genome = 'hg19', name = 'TARID', type='p',frame=TRUE,col.frame='black')

tracks[['TARID']]=eqtl_track


pdf(sprintf('%s/tarid_tcf21_coloc.pdf',fig_dir,chr,start,end),height=15,width=8)
plotTracks(tracks[c(1:11,14,15,12,13)], chromosome = chr, from = start, to = end)
dev.off()