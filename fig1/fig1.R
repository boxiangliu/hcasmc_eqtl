library(data.table)
library(ggrepel)
library(cowplot)
library(ggplot2)
library(grid)
library(gridExtra)
library(VennDiagram)
library(stringr)
library(Gviz)
library(rtracklayer)
library(gtable)

out_dir='../figures/fig1/'
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# Fig. 1A:
load('../processed_data/160603/mds.Rdata')
fig1a=p3


# Fig. 1B:
temp=readRDS('../processed_data/atacseq_similarity/jaccard_similarity/jaccard.rds')
jaccard=temp[[1]]
color_map=temp[[2]]


show_n=5
jaccard[,label:=ifelse( (rank<=show_n) | (rank>(length(jaccard)-show_n)),as.character(sample),'')]
# jaccard[rank<70,list(sample,gtex)]
top5_label=jaccard[rank<=show_n,as.character(sample)]
bottom5_label=jaccard[rank>=(length(jaccard)-show_n),as.character(sample)]

fig1b=ggplot(jaccard,aes(sample,jaccard,label=label,color=gtex))+
	geom_point(size=1)+
	theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
	xlab(sprintf('ENCODE Sample (n=%i)',nrow(jaccard)))+
	ylab('Epigenomic Similarity\n(Jaccard Index)')+
	annotate('text',x=20,y=0.10,label=paste(c('Top 5',top5_label),collapse='\n'))+
	annotate('text',x=80,y=0.19,label=paste(c('Bottom 5',bottom5_label),collapse='\n'))+
	scale_color_manual(values=color_map,guide='none')

# Fig. 1C:
esi=fread('../processed_data/160715/esi.quant_norm.all_tissues.txt')

annotation=fread('/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.bed',select=c(5,6,7),col.names=c('gene_id','gene_name','type'))


esi[,list(gene,HCASMC)]
esi_pc=esi[gene %in% annotation[type=='protein_coding',gene_id]]
esi_hcasmc=esi_pc[,list(gene,HCASMC)]
setorder(esi_hcasmc,-HCASMC)

temp=readRDS(sprintf('../processed_data/160715/find_hcasmc_specific_genes.info_theory.quant_norm/rpkm.rds'))
rpkmn=temp[[2]]
tissue_median=temp[[3]]
rowdata=temp[[4]]
coldata=temp[[5]]


tissue_median$gene_id=rowdata$Name
tissue_median$gene_name=rowdata$Description
stopifnot(all(rownames(tissue_median)==tissue_median$gene_id))

idx=match(esi_hcasmc[1:15,gene],tissue_median$gene_id)
hcasmc_spec_gene=tissue_median[idx,]


hcasmc=hcasmc_spec_gene[,c('gene_id','gene_name','HCASMC')]


gtex=hcasmc_spec_gene[,-which(colnames(hcasmc_spec_gene)=='HCASMC')]
gtex_long=melt(gtex,id.vars=c('gene_id','gene_name'),variable.name='tissue',value.name='rpkm')

x_axis=hcasmc_spec_gene$gene_name[order(-hcasmc_spec_gene$HCASMC)]

fig1c=ggplot(gtex_long,aes(x=gene_name,y=rpkm))+
	geom_violin()+
	scale_y_log10()+
	geom_point(data=hcasmc,aes(x=gene_name,y=HCASMC),color='#FF0066')+
	lims(x=x_axis)+coord_flip()+
	theme(axis.title.y=element_blank(),axis.text.y=element_text(size=8))

# Fig. 1D:
df1=data.frame(
	Category=c(
		'Cell division and\nproliferation',
		'Phenotype\ntransition',
		'ECM\nsecretion'),
	`Upregulated Pathways`=c(
		'G2M checkpoint, mitotic spindle\nMYC target, E2F targets\nmTORC1, DNA repair',
		'Epithelial mesencymal transition',
		'Protein secretion\nUnfolded protein response'))
g=tableGrob(df1, rows = NULL)
g=gtable_add_grob(g,grobs=rectGrob(gp=gpar(fill=NA,lwd=2)),t=2,b=nrow(g),l=1,r=ncol(g))
g=gtable_add_grob(g,grobs=rectGrob(gp=gpar(fill=NA,lwd=2)),t=1,l=1,r=ncol(g))
fig1d=gtable_add_grob(g,grobs=rectGrob(gp=gpar(fill=NA,lwd=2)),t=1,b=nrow(g),l=2)



# Fig. 1E:
intersect=fread('../processed_data/hcasmc_specific_open_chromatin/intersect/intersect.bed')
hcasmc_col=which(names(intersect)=='HCASMC')
idx=which(intersect[,hcasmc_col,with=FALSE]==1)
n_hcasmc_peaks=length(idx)
n_hcasmc_spec_peaks=sum(intersect[idx,4,with=FALSE]==1)

fig1e=draw.pairwise.venn(
	area1=n_hcasmc_peaks, 
	area2=nrow(intersect),
	cross.area=n_hcasmc_peaks-n_hcasmc_spec_peaks,
	category=c("HCASMC", "ENCODE"),
	lty=rep("blank", 2), 
	fill=c("blue", "red"), 
	alpha = rep(0.5, 2), 
	cat.pos = c(0, 0), 
	cat.dist = rep(0.025, 2),
	cex=1.3,
	cat.cex=1.5)

# Fig 1F:
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

combine_data=function(bigwig_fn,tissue,chr,start,end){
	
	chr_=chr
	start_=start 
	end_=end

	container=list()

	for (i in seq(length(bigwig_fn))){
		sample=tissue[i]
		f=bigwig_fn[i]
		print(sample)
		print(f)
		bw=get_bw(f,seq=chr_,start_,end_)
		bw$rsid='.'
		bw$data=sample
		container[[f]]=bw
	}
	bw=Reduce(rbind,container)

	to_plot=rbind(bw[,list(chr=seqnames,pos,logp=score,rsid,data)])

	return(to_plot)
}

metadata=fread('../data/encode/dnase_seq_bigwig/metadata.tsv',
	select=c(1,3,4,7,8,38),
	col.names=c('file_accession',
				'output_type',
				'experiment_accession',
				'biosample_term_name',
				'biosample_type',
				'size'))

in_dir='../data/encode/dnase_seq_bigwig/'

metadata[biosample_term_name %in% 
	c('fibroblast of lung',
	  'bronchial epithelial cell',
	  'astrocyte',
	  'skeletal muscle cell',
	  'stromal cell of bone marrow',
	  'keratinocyte',
	  'B cell',
	  'iris pigment epithelial cell') & 
	output_type=='raw signal' &
	size==max_size,biosample_term_name,by='biosample_term_name']

metadata=metadata[output_type == 'raw signal' & biosample_term_name %in% 
	c('fibroblast of lung',
	  'bronchial epithelial cell',
	  'astrocyte',
	  'skeletal muscle cell',
	  'stromal cell of bone marrow',
	  'keratinocyte',
	  'B cell',
	  'iris pigment epithelial cell')][,max_size:=max(size),by='biosample_term_name']
metadata=metadata[size==max_size]


bigwig_fn=metadata[output_type=='raw signal',sprintf('%s/%s.bigWig',in_dir,file_accession)]
bigwig_fn=c('../data/atacseq/fbs/2305/out/signal/macs2/pooled_rep/CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.nodup.tn5_pooled.pf.fc.signal.bigwig',
			bigwig_fn)
tissue=metadata[output_type=='raw signal',biosample_term_name]
tissue=c('HCASMC',tissue)

chr='chr10'
start=10838e4
end=10843e4
to_plot=unique(combine_data(bigwig_fn,tissue,chr,start,end))


fig1f=ggplot(to_plot,aes(pos,logp))+facet_grid(data~.,scales='free_y')+
	geom_area()+theme_bw()+
	geom_text(aes(label=data),x=Inf,y=Inf,hjust=1.0,vjust=1.0)+
	xlab(chr)+ylab('Signal')+
	theme(axis.title=element_blank(),
		axis.text=element_blank(),
		panel.grid=element_blank(),
		panel.border=element_blank(),
		axis.ticks=element_blank(),
		panel.spacing.y=unit(0,'cm'),
		strip.background=element_blank(),
		strip.text.y=element_blank())




# Fig. 1:
pdf(sprintf('%s/fig1.pdf',out_dir),height=8,width=16)
grid.arrange(arrangeGrob(fig1a,
	arrangeGrob(fig1b,
		arrangeGrob(fig1c,fig1d,ncol=2),
		arrangeGrob(gTree(children=fig1e),fig1f,ncol=2),
		ncol=1),ncol=2))
dev.off()