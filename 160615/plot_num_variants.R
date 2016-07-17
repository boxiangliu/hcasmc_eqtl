# command args: 
args=commandArgs(T,T)
ld_file=args$ld
bim_file=args$bim
variants_file=args$variants
ld_thresh=as.numeric(args$ld)
figure_file=args$figure

# ld_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/LD_fdr0.1_dprime0.8/target_variants.ld'
# bim_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/LD_fdr0.1_dprime0.8/target_variants.bim'
# variants_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/LD_fdr0.1_dprime0.8/combined.txt'
# figure_file='/srv/persistent/bliu2/HCASMC_eQTL/figures/160615/num_ld_variant_by_tissue.dprime0.8.fdr0.1.ld0.2.pdf'
# ld_thresh=0.2

# read ld file
ld=fread(ld_file,header=F)
bim=fread(bim_file)
variants=fread(variants_file,header=F)


# create id from chr and pos:
bim[,id:=paste(V1,V4,sep='_')]


# cast ld into a matrix:
ld=as.matrix(ld)
colnames(ld)=bim$id
rownames(ld)=bim$id


# create id from chr and pos: 
setnames(variants,c('chr','pos','gene_id','variant_id','type'))
variants[,id:=paste(chr,pos,sep='_')]


# subset LD matrix to GWAS columns: 
gwas=ld[,colnames(ld)%in%variants[type=='GWAS',id]]


# calcualte max LD per row/eQTL: 
max_ld=c()
for (i in 1:nrow(gwas)){
	max_ld=c(max_ld,max(gwas[i,],na.rm=T))
}
gwas_max_ld=data.table(id=rownames(gwas),max_ld=max_ld)


# high LD variants for HCASMC: 
gwas_max_ld[id%in%variants[type=='HCASMC',id]&max_ld>0.8,]

# Two high LD variants: 
#              id   max_ld
# 1: 10_104611764 0.985949
# 2:   17_2208874 1.000000

# Corresponding eGene:
# 10_104611764 ENSG00000214435.3 AS3MT
# 17_2208874 ENSG00000167720.8 SRR

# Gene function: 
# AS3MT catalyzes the transfer of a methyl group from S-adenosyl-L-methionine (AdoMet) to trivalent arsenical and may play a role in arsenic metabolism (Lin et al., 2002 [PubMed 11790780]).[supplied by OMIM, Mar 2008]
# SRR (Serine Racemase) generates D-serine from L-serine. Diseases associated with SRR include streptococcal meningitis and serine deficiency. Among its related pathways are Glycine, serine and threonine metabolism and serine and glycine biosynthesis. GO annotations related to this gene include calcium ion binding and magnesium ion binding.


# subset for variants with LD > 0.2:
tissue_max_ld=data.table()
for (t in unique(variants$type)){
	if (t=='GWAS') next()
	idx=with(gwas_max_ld, (id%in%variants[type==t,id] & max_ld>ld_thresh))
	temp=data.table(type=t,gwas_max_ld[idx,])
	tissue_max_ld=rbind(tissue_max_ld,temp)
}
tissue_max_ld[,num:=.N,by=type]


# plot the number of high LD eQTLs:
tissue_max_ld$type=as.factor(tissue_max_ld$type)
p=ggplot(tissue_max_ld,aes(x=reorder(type,-num,FUN=mean),y=num))+geom_point()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+xlab('Tissue')+ylab('Number of variants')+background_grid(major='x')
save_plot(figure_file,p,base_height=7, base_width=7)
