library(data.table)
library(cowplot)

# Variables: 
fig_dir='../figures/hcasmc_specific_open_chromatin/classify_peaks/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir)}

# Read gencode: 
gencode=fread('../../shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf')


# Give meaningful names and remove uncessary columns:
setnames(gencode,c('V1','V3','V4','V5','V7'),c('chr','type','start','end','strand'))
gencode[,c('V2','V6','V8','V9'):=NULL]


# Subset to genes:
gencode=gencode[type=='gene']


# Get the TSS: 
gencode[,tss:=ifelse(strand=='+',start,end)]


# Get HCASMC peaks (raw peaks, not summit extend peaks):
peaks=fread('zcat ../processed_data/hcasmc_specific_open_chromatin/encode_plus_hcasmc/2305_ppr.IDR0.1.filt.narrowPeak.gz')


# Give meaningfulnames and remove uncessary columms:
setnames(peaks,c('V1','V2','V3','V10'),c('chr','start','end','summit'))
peaks[,c('V4','V5','V6','V7','V8','V9'):=NULL]


# Keep only autosomes: 
peaks=peaks[(chr!='chrX')&(chr!='chrY')&(chr!='chrM'),]


# Calculate peak size:
peaks[,size:=end-start]


# Calculate summit extend summit boundary:
peaks[,c('start_se','end_se'):=list(start+summit-75,start+summit+75)]


# Assign each peak a unique ID:
peaks[,id:=paste(chr,start,end,summit,sep='_')]


# Overlap HCASMC peaks with gencode TSS:
gencode[,c('tss_start','tss_end'):=tss]
setkey(peaks,chr,start,end)
setkey(gencode,chr,tss_start,tss_end)
tss_ol=foverlaps(gencode,peaks,nomatch=0)
tss_ol=unique(tss_ol[,.(chr,start,end,summit,size,start_se,end_se,id)])


# Overlap HCASMC peaks with gencode gene body:
setkey(gencode,chr,start,end)
gene_body_ol=foverlaps(peaks,gencode,type='any',nomatch=0)
gene_body_ol=unique(gene_body_ol[,.(chr,i.start,i.end,summit,size,start_se,end_se,id)])
setnames(gene_body_ol,c('i.start','i.end'),c('start','end'))
gene_body_ol=gene_body_ol[!id%in%tss_ol$id,]


# Get intergenic peaks:
intergenic_ol=peaks[(!id%in%gene_body_ol$id)&(!id%in%tss_ol$id),]


# Plot peak size vs peak type:
to_plot=rbind(data.table(size=unique(intergenic_ol[,.(chr,start,end,size)])$size,type='intergenic'),
			  data.table(size=unique(tss_ol[,.(chr,start,end,size)])$size,type='tss'),
			  data.table(size=unique(gene_body_ol[,.(chr,start,end,size)])$size,type='gene_body'),
			  data.table(size=unique(peaks[,.(chr,start,end,size)])$size,type='all'))
to_plot$type=factor(to_plot$type,levels=c('all','tss','gene_body','intergenic'))
p1=ggplot(to_plot,aes(x=type,y=size))+geom_boxplot(outlier.size=-1)+ylim(0,2.3e3)
save_plot(sprintf('%s/peak_size_vs_type.pdf',fig_dir),p1)


# Get HCASMC specific peaks:
peak_specific=fread('../processed_data/hcasmc_specific_open_chromatin/peak_specificity_filt/HCASMC.bed')


# 
tss_ol_specific=merge(tss_ol,peak_specific,by.x=c('chr','start_se','end_se'),by.y=c('chr','start','end'))
tss_ol_specific=tss_ol_specific[,list(psi=min(psi)),by=c('chr','start','end')]


gene_body_ol_specific=merge(gene_body_ol,peak_specific,by.x=c('chr','start_se','end_se'),by.y=c('chr','start','end'))
gene_body_ol_specific=gene_body_ol_specific[,list(psi=min(psi)),by=c('chr','start','end')]


intergenic_ol_specific=merge(intergenic_ol,peak_specific,by.x=c('chr','start_se','end_se'),by.y=c('chr','start','end'))
intergenic_ol_specific=intergenic_ol_specific[,list(psi=min(psi)),by=c('chr','start','end')]

peak_ol_specific=merge(peaks,peak_specific,by.x=c('chr','start_se','end_se'),by.y=c('chr','start','end'))
peak_ol_specific=peak_ol_specific[,list(psi=min(psi)),by=c('chr','start','end')]

to_plot=rbind(data.table(psi=intergenic_ol_specific$psi,type='intergenic'),
			  data.table(psi=tss_ol_specific$psi,type='tss'),
			  data.table(psi=gene_body_ol_specific$psi,type='gene_body'),
			  data.table(psi=peak_ol_specific$psi,type='all'))
to_plot$type=factor(to_plot$type,levels=c('all','tss','gene_body','intergenic'))
p2=ggplot(to_plot,aes(type,psi))+geom_violin()
save_plot(sprintf("%s/peak_specificity_vs_type.pdf",fig_dir),p2)
