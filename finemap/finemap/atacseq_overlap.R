library(data.table)
library(stringr)
library(cowplot)
library(rtracklayer)
library(ggrepel)

atac_dir='../data/atacseq/fbs/'
bw_fn='../data/atacseq/fbs/2305/out/signal/macs2/pooled_rep/CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.nodup.tn5_pooled.pf.fc.signal.bigwig'
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

combine_data=function(bigwig_fn,nikpay,howson,eqtl,chr,start,end){
	
	chr_=chr
	start_=start 
	end_=end


	bw=get_bw(bigwig_fn,seq=chr_,start_,end_)
	bw$rsid='.'

	ni=nikpay[chr==chr_&pos>=start_&pos<=end_,]

	ho=howson[chr==chr_&pos>=start_&pos<=end_,]

	eq=eqtl[chr==chr_&pos>=start_&pos<=end_,]

	to_plot=rbind(bw[,list(chr=seqnames,pos,logp=score,rsid,data='ATACseq')],
		ni[,list(chr,pos,logp,rsid,data='Nikpay')],
		ho[,list(chr,pos,logp,rsid,data='Howson')],
		eq[,list(chr,pos,logp,rsid,data='RASQUAL')])

	return(to_plot)
}

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


# TCF21:
chr='chr6'
start=134e6
end=134.4e6
eqtl_fn='../processed_data/rasqual/output_pval/chr6/ENSG00000118526.6_TCF21.pval.txt'

to_plot=combine_data(bw_fn,nikpay,howson,eqtl_fn,chr,start,end)

p1=ggplot(to_plot,aes(pos,logp))+facet_grid(data~.,scales='free_y')+
	geom_point(data=to_plot[data%in%c('Nikpay','Howson','RASQUAL'),])+
	geom_area(data=to_plot[data=='ATACseq',])+theme_bw()

to_plot2=to_plot[pos>=134209800&pos<=134215600]
p2=ggplot(to_plot2,aes(pos,logp))+
	facet_grid(data~.,scales='free_y')+
	geom_point(data=to_plot2[data%in%c('Nikpay','Howson','RASQUAL'),])+
	geom_area(data=to_plot2[data=='ATACseq',])+theme_bw()+
	geom_text(aes(label=ifelse(pos%in%c(134209837,134214525),rsid,'')),hjust=-0.1)

p3=plot_ase(eq[rsid=='rs2327429',],'rs2327429')
p4=plot_ase(eq[rsid=='rs12190287',],'rs2327429')

pdf(sprintf('%s/tcf21.pdf',fig_dir))
p1;p2;p3;p4
dev.off()


# FES:
chr='chr15'
start=905e5
end=920e5
eqtl_fn='../processed_data/rasqual/output_pval/chr15/ENSG00000182511.7_FES.pval.txt'
eqtl=fread(eqtl_fn)
eqtl[,logp:=-log10(pval)]


to_plot=combine_data(bw_fn,nikpay,howson,eqtl,chr,start,end)

subset=to_plot
p1=ggplot(subset,aes(pos,logp))+facet_grid(data~.,scales='free_y')+
	geom_point(data=subset[data%in%c('Nikpay','Howson','RASQUAL'),],alpha=0.5)+
	geom_area(data=subset[data=='ATACseq',])+theme_bw()

subset=to_plot[pos>=913e5&pos<=915e5,]
p2=ggplot(subset,aes(pos,logp))+facet_grid(data~.,scales='free_y')+
	geom_point(data=subset[data%in%c('Nikpay','Howson','RASQUAL'),],alpha=0.5)+
	geom_area(data=subset[data=='ATACseq',])+theme_bw()


subset=to_plot[pos>=91426300&pos<=91448950,]
p3=ggplot(subset,aes(pos,logp,label=ifelse(rsid=='rs2521501',rsid,'')))+facet_grid(data~.,scales='free_y')+
	geom_point(data=subset[data%in%c('Nikpay','Howson','RASQUAL'),],alpha=0.5)+
	geom_area(data=subset[data=='ATACseq',])+theme_bw()+geom_text(hjust=1.1)


snps=unlist(subset[data=='RASQUAL'&logp>=7,rsid])
p4=ggplot(subset,aes(pos,logp,label=ifelse(rsid%in%snps,rsid,'')))+facet_grid(data~.,scales='free_y')+
	geom_point(data=subset[data%in%c('Nikpay','Howson','RASQUAL'),],alpha=0.5)+
	geom_area(data=subset[data=='ATACseq',])+geom_text_repel(force=3)+theme_bw()

container=list()
for (snp in snps){
	p=plot_ase(eqtl[rsid==snp,],snp)
	container[[snp]]=p
}

pdf(sprintf('%s/fes.pdf',fig_dir))
p1;p2;p3;p4
for (i in seq_along(container)){print(container[[i]])}
dev.off()


# SIPA1:
chr='chr11'
start=65e6
end=66e6
eqtl_fn='../processed_data/rasqual/output_pval/chr11/ENSG00000213445.4_SIPA1.pval.txt'
eqtl=fread(eqtl_fn)
eqtl[,logp:=-log10(pval)]

to_plot=combine_data(bw_fn,nikpay,howson,eqtl,chr,start,end)

subset=to_plot
p1=ggplot(subset,aes(pos,logp))+facet_grid(data~.,scales='free_y')+
	geom_point(data=subset[data%in%c('Nikpay','Howson','RASQUAL'),],alpha=0.2)+
	geom_area(data=subset[data=='ATACseq',])+theme_bw()


subset=to_plot[pos>=65375e3&pos<=65475e3]
snps=eqtl[logp==max(logp),rsid]
snps=c(snps,eqtl[logp>=6&pos<=65412500&pos>=65400000,rsid])
p2=ggplot(subset,aes(pos,logp,label=ifelse(rsid%in%snps,rsid,'')))+facet_grid(data~.,scales='free_y')+
	geom_point(data=subset[data%in%c('Nikpay','Howson','RASQUAL'),],alpha=0.2)+
	geom_area(data=subset[data=='ATACseq',])+theme_bw()+geom_text_repel(force=3)


subset=to_plot[pos>=65400000&pos<=65412500]
snps=eqtl[logp>=6&pos<=65412500&pos>=65400000,rsid]
p3=ggplot(subset,aes(pos,logp,label=ifelse(rsid%in%snps,rsid,'')))+facet_grid(data~.,scales='free_y')+
	geom_point(data=subset[data%in%c('Nikpay','Howson','RASQUAL'),],alpha=0.2)+
	geom_area(data=subset[data=='ATACseq',])+theme_bw()+geom_text_repel(force=3)


subset=to_plot[pos>=65380000&pos<=65400000]
snps=eqtl[logp==max(logp),rsid]
p4=ggplot(subset,aes(pos,logp,label=ifelse(rsid%in%snps,rsid,'')))+facet_grid(data~.,scales='free_y')+
	geom_point(data=subset[data%in%c('Nikpay','Howson','RASQUAL'),],alpha=0.2)+
	geom_area(data=subset[data=='ATACseq',])+theme_bw()+geom_text_repel(force=3)


snps=eqtl[logp==max(logp),rsid]
snps=c(snps,eqtl[logp>=6&pos<=65412500&pos>=65400000,rsid])

container=list()
for (snp in snps){
	p=plot_ase(eqtl[rsid==snp,],snp)
	container[[snp]]=p
}


pdf(sprintf('%s/sipa1.pdf',fig_dir))
p1;p2;p3;p4
for (i in seq_along(container)){print(container[[i]])}
dev.off()


# MRAS:
chr='chr3'
start=1.375e8
end=1.385e8
eqtl_fn='../processed_data/rasqual/output_pval/chr3/ENSG00000158186.8_MRAS.pval.txt'
eqtl=fread(eqtl_fn)
eqtl[,logp:=-log10(pval)]

to_plot=combine_data(bw_fn,nikpay,howson,eqtl,chr,start,end)

subset=to_plot
p1=ggplot(subset,aes(pos,logp))+facet_grid(data~.,scales='free_y')+
	geom_point(data=subset[data%in%c('Nikpay','Howson','RASQUAL'),],alpha=0.2)+
	geom_area(data=subset[data=='ATACseq',])+theme_bw()

subset_cast=dcast(subset[data!='ATACseq'],chr+pos~data,value.var='logp')
p2=ggplot(subset_cast,aes(Nikpay,RASQUAL))+geom_point(size=2,alpha=0.5)
p3=ggplot(subset_cast,aes(Howson,RASQUAL))+geom_point(size=2,alpha=0.5)

subset=to_plot[pos>=138000000&pos<=138250000]

snps=c(subset[data=='Nikpay'&logp>7.5,pos],
	subset[data=='Howson'&logp>6.8,pos])

p4=ggplot(subset,aes(pos,logp,label=ifelse(pos%in%snps,rsid,'')))+facet_grid(data~.,scales='free_y')+
	geom_point(data=subset[data%in%c('Nikpay','Howson','RASQUAL'),],alpha=0.2)+
	geom_area(data=subset[data=='ATACseq',])+theme_bw()+geom_text_repel()

subset=to_plot[pos>=138075000&pos<=138100000,]
p5=ggplot(subset,aes(pos,logp,label=ifelse(pos%in%snps,rsid,'')))+facet_grid(data~.,scales='free_y')+
	geom_point(data=subset[data%in%c('Nikpay','Howson','RASQUAL'),],alpha=0.2)+
	geom_area(data=subset[data=='ATACseq',])+theme_bw()+geom_text_repel()


subset=to_plot[pos>=138100000&pos<=138125000,]
p6=ggplot(subset,aes(pos,logp,label=ifelse(pos%in%snps,rsid,'')))+facet_grid(data~.,scales='free_y')+
	geom_point(data=subset[data%in%c('Nikpay','Howson','RASQUAL'),],alpha=0.2)+
	geom_area(data=subset[data=='ATACseq',])+theme_bw()+geom_text_repel()


container=list()
for (snp in snps){
	rsid=eqtl[pos==snp,rsid]
	if (length(rsid)!=0){
		p=plot_ase(eqtl[pos==snp,],rsid)
		container[[rsid]]=p
	}
}

pdf(sprintf('%s/mras.pdf',fig_dir))
p1;p2;p3;p4;p5;p6
for (i in seq_along(container)) {print(container[[i]])}
dev.off()


# AS3MT:
chr='chr10'
start=1.043e8
end=1.053e8
eqtl_fn='../processed_data/rasqual/output_pval/chr10/ENSG00000214435.3_AS3M.pval.txt'
eqtl=fread(eqtl_fn)
eqtl[,logp:=-log10(pval)]

to_plot=combine_data(bw_fn,nikpay,howson,eqtl,chr,start,end)
to_plot=unique(to_plot)

subset=to_plot
p1=ggplot(subset,aes(pos,logp))+facet_grid(data~.,scales='free_y')+
	geom_point(data=subset[data%in%c('Nikpay','Howson','RASQUAL'),],alpha=0.2)+
	geom_area(data=subset[data=='ATACseq',])+theme_bw()

subset=to_plot[pos<=104750000&pos>=104500000]
p2=ggplot(subset,aes(pos,logp))+facet_grid(data~.,scales='free_y')+
	geom_point(data=subset[data%in%c('Nikpay','Howson','RASQUAL'),],alpha=0.2)+
	geom_area(data=subset[data=='ATACseq',])+theme_bw()


subset_cast=dcast(subset[data!='ATACseq'],chr+pos+rsid~data,value.var='logp')
p3=ggplot(subset_cast,aes(Nikpay,RASQUAL))+geom_point(size=2,alpha=0.5)+theme_bw()
p4=ggplot(subset_cast,aes(Howson,RASQUAL))+geom_point(size=2,alpha=0.5)+theme_bw()


pdf(sprintf('%s/as3mt.pdf',fig_dir))
p1;p2;p3;p4
dev.off()


# SMAD3:
chr='chr15'
start=67e6
end=67.5e6
eqtl_fn='../processed_data/rasqual/output_pval/chr15/ENSG00000166949.11_SMAD3.pval.txt'
eqtl=fread(eqtl_fn)
eqtl[,logp:=-log10(pval)]

to_plot=combine_data(bw_fn,nikpay,howson,eqtl,chr,start,end)
to_plot=unique(to_plot)

subset=to_plot
subset_cast=dcast(subset[data!='ATACseq'],chr+pos+rsid~data,value.var='logp',fun.aggregate=max)
snps=subset_cast[Nikpay>=7.5&RASQUAL>=10,rsid]

p1=ggplot(subset,aes(pos,logp,label=ifelse(rsid%in%snps,rsid,'')))+facet_grid(data~.,scales='free_y')+
	geom_point(data=subset[data%in%c('Nikpay','Howson','RASQUAL'),],alpha=0.2)+
	geom_area(data=subset[data=='ATACseq',])+theme_bw()+geom_text_repel()


subset=to_plot[pos>=67440750&pos<=67456630]
p2=ggplot(subset,aes(pos,logp,label=ifelse(rsid%in%snps,rsid,'')))+facet_grid(data~.,scales='free_y')+
	geom_point(data=subset[data%in%c('Nikpay','Howson','RASQUAL'),],alpha=0.5)+
	geom_area(data=subset[data=='ATACseq',])+theme_bw()+geom_text_repel()


container=list()
for (snp in snps){
	p=plot_ase(eqtl[rsid==snp,],snp)
	container[[snp]]=p
}

pdf(sprintf('%s/smad3.pdf',fig_dir))
p1;p2
for (i in seq_along(container)){print(container[[i]])}
dev.off()


# RP11-54A9.1:
chr='chr12'
start=7.6e7
end=7.7e7
eqtl_fn='../processed_data/rasqual/output_pval/chr12/ENSG00000257219.1_RP11-54A9.1.pval.txt'
eqtl=fread(eqtl_fn)
eqtl[,logp:=-log10(pval)]

to_plot=combine_data(bw_fn,nikpay,howson,eqtl,chr,start,end)
to_plot=unique(to_plot)

subset=to_plot

subset_cast=dcast(subset[data!='ATACseq'],chr+pos+rsid~data,value.var='logp',fun.aggregate=max)
snps=subset_cast[Nikpay>=4.5|RASQUAL>=7,rsid]


p1=ggplot(subset,aes(pos,logp,label=ifelse(rsid%in%snps,rsid,'')))+facet_grid(data~.,scales='free_y')+
	geom_point(data=subset[data%in%c('Nikpay','Howson','RASQUAL'),],alpha=0.2)+
	geom_area(data=subset[data=='ATACseq',])+theme_bw()+geom_text_repel()

p2=ggplot(subset_cast,aes(Nikpay,RASQUAL,label=ifelse(rsid%in%snps,rsid,'')))+geom_point(size=2,alpha=0.5)+theme_bw()+geom_text_repel()


pdf(sprintf('%s/rp11-54a9.1.pdf',fig_dir))
p1;p2
dev.off()


# PDLIM5:
chr='chr4'
start=95e6
end=96e6
eqtl_fn='../processed_data/rasqual/output_pval/chr4/ENSG00000163110.10_PDLIM5.pval.txt'
eqtl=fread(eqtl_fn)
eqtl[,logp:=-log10(pval)]

to_plot=combine_data(bw_fn,nikpay,howson,eqtl,chr,start,end)
to_plot=unique(to_plot)

subset=to_plot
p1=ggplot(subset,aes(pos,logp,label=ifelse(rsid%in%snps,rsid,'')))+facet_grid(data~.,scales='free_y')+
	geom_point(data=subset[data%in%c('Nikpay','Howson','RASQUAL'),],alpha=0.2)+
	geom_area(data=subset[data=='ATACseq',])+theme_bw()+geom_text_repel()

subset_cast=dcast(subset[data!='ATACseq'],chr+pos+rsid~data,value.var='logp',fun.aggregate=max)
p2=ggplot(subset_cast,aes(Nikpay,RASQUAL,label=ifelse(rsid%in%snps,rsid,'')))+geom_point(size=2,alpha=0.5)+theme_bw()+geom_text_repel()


pdf(sprintf('%s/pdlim5.pdf',fig_dir))
p1;p2
dev.off()