library(data.table)
library(cowplot)
library(limma)

fig_dir='../figures/compare_rasqual_and_fastqtl/gtex_replication/'
if (!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}


# eGene replication: 
ras_fn='../processed_data/rasqual/output_merged/treeQTL/eGenes.tsv'
fas_fn='../processed_data/eqtl/fastqtl/output/treeQTL/eGenes.tsv'
gtex_fn='../processed_data/compare_rasqual_and_fastqtl/gtex_eqtl/TreeQTL/eGenes.tsv'

ras=fread(ras_fn)
fas=fread(fas_fn)
gtex=fread(gtex_fn)


y=data.frame()
for (i in c(seq(10,160,10),seq(260,1220,100))){
	ras_egenes=ras[1:i,family]
	fas_egenes=fas[1:i,family]
	ras_rep_pval=-log10(gtex[family%in%ras_egenes,fam_p])
	fas_rep_pval=-log10(gtex[family%in%fas_egenes,fam_p])
	x=rbind(data.frame(pval=ras_rep_pval,method='ras',n=i,n_rep=length(ras_rep_pval)),data.frame(pval=fas_rep_pval,method='fas',n=i,n_rep=length(fas_rep_pval)))
	y=rbind(y,x)
}
y$n=as.factor(y$n)
sum(ras$family%in%gtex$family) # 1131
sum(fas$family%in%gtex$family) # 162


p1=ggplot(y,aes(n,pval,fill=method))+geom_boxplot()+ggtitle('eGenes')+theme(axis.text.x=element_text(angle=90))
p1_2=ggplot(y,aes(n,n_rep,color=method))+geom_point()+ggtitle('eGenes')+theme(axis.text.x=element_text(angle=90))+scale_y_continuous(breaks=seq(0,1200,100))
pdf(sprintf('%s/eGenes.pdf',fig_dir));p1;p1_2;dev.off()


merged=merge(fas[,list(family,fam_p)],ras[,list(family,fam_p)],,by='family',all=TRUE,suffixes = c(".fas",".ras"))
c=merged[,list(ras=!is.na(fam_p.ras),fas=!is.na(fam_p.fas))]
pdf(sprintf('%s/eGenes.venn.pdf',fig_dir))
a=vennCounts(c)
vennDiagram(a)
a[,'Counts']=round(a[,'Counts']/sum(a[,'Counts']),4)
vennDiagram(a)
dev.off()


# eSNP replication:
ras_fn='../processed_data/rasqual/output_merged/treeQTL/eSNPs.tsv'
fas_fn='../processed_data/eqtl/fastqtl/output/treeQTL/eSNPs.tsv'
gtex_fn='../processed_data/compare_rasqual_and_fastqtl/gtex_eqtl/TreeQTL/eSNPs.tsv'

ras=fread(ras_fn)
fas=fread(fas_fn)
gtex=fread(gtex_fn)

y=data.frame()
for (i in c(seq(1000,13000,1000),seq(23000,210000,10000))){
	ras_eSNPs=ras[1:i,family]
	fas_eSNPs=fas[1:min(i,nrow(fas)),family]
	ras_rep_pval=-log10(gtex[family%in%ras_eSNPs,fam_p])
	fas_rep_pval=-log10(gtex[family%in%fas_eSNPs,fam_p])
	x=rbind(data.frame(pval=ras_rep_pval,method='ras',n=i,n_rep=length(ras_rep_pval)),data.frame(pval=fas_rep_pval,method='fas',n=i,n_rep=length(fas_rep_pval)))
	y=rbind(y,x)
}
y$n=as.factor(y$n)
sum(ras$family%in%gtex$family) # 1131
sum(fas$family%in%gtex$family) # 162


p2=ggplot(y,aes(n,pval,fill=method))+geom_boxplot()+ggtitle('eSNPs')+theme(axis.text.x=element_text(angle=90))
p2_2=ggplot(y,aes(n,n_rep,color=method))+geom_point()+ggtitle('eSNPs')+theme(axis.text.x=element_text(angle=90))
pdf(sprintf('%s/eSNPs.pdf',fig_dir));p2;p2_2;dev.off()


merged=merge(fas[,list(family,fam_p)],ras[,list(family,fam_p)],,by='family',all=TRUE,suffixes = c(".fas",".ras"))
c=merged[,list(ras=!is.na(fam_p.ras),fas=!is.na(fam_p.fas))]
pdf(sprintf('%s/eSNPs.venn.pdf',fig_dir))
a=vennCounts(c)
vennDiagram(a)
a[,'Counts']=round(a[,'Counts']/sum(a[,'Counts']),4)
vennDiagram(a)
dev.off()

# eAssociation replication: 
ras_fn='../processed_data/rasqual/output_merged/treeQTL/eAssoc_eSNPs.tsv'
fas_fn='../processed_data/eqtl/fastqtl/output/treeQTL/eAssoc_eSNPs.tsv'
gtex_fn='../processed_data/compare_rasqual_and_fastqtl/gtex_eqtl/TreeQTL/eAssoc_eSNPs.tsv'

ras=fread(ras_fn)
fas=fread(fas_fn)
gtex=fread(gtex_fn)


# Add ID column: 
gtex[,id:=paste(gene,SNP,sep='_')]
ras[,id:=paste(gene,SNP,sep='_')]
fas[,id:=paste(gene,SNP,sep='_')]

y=data.frame()
for (i in c(seq(2000,20000,2000),seq(40000,240000,20000))){
	ras_eAssoc=ras[1:i,id]
	fas_eAssoc=fas[1:min(i,nrow(fas)),id]
	ras_rep_pval=-log10(gtex[id%in%ras_eAssoc,BBFDR])
	fas_rep_pval=-log10(gtex[id%in%fas_eAssoc,BBFDR])
	x=rbind(data.frame(pval=ras_rep_pval,method='ras',n=i,n_rep=length(ras_rep_pval)),data.frame(pval=fas_rep_pval,method='fas',n=i,n_rep=length(fas_rep_pval)))
	y=rbind(y,x)
}
y$n=as.factor(y$n)
sum(ras$id%in%gtex$id) # 94264
length(ras$id) # 248421
sum(fas$id%in%gtex$id) # 18253
length(fas$id) # 20676


p3=ggplot(y,aes(n,pval,fill=method))+geom_boxplot()+ggtitle('eAssociations')+theme(axis.text.x=element_text(angle=90))
p3_2=ggplot(y,aes(n,n_rep,color=method))+geom_point()+ggtitle('eAssociations')+theme(axis.text.x=element_text(angle=90))
pdf(sprintf('%s/eAssociations.pdf',fig_dir));p3;p3_2;dev.off()


# Generate sensitivity vs FDR data: 
FDR=seq(0,0.05,1e-4)
sens=matrix(NA,nrow=length(FDR),ncol=2)
colnames(sens)=c('ras','fas')


for (i in 1:length(FDR)){
	sub=ras[BBFDR<FDR[i],id]
	sens[i,'ras']=length(intersect(gtex$id,sub))# takes a while

	sub=fas[BBFDR<FDR[i],id]
	sens[i,'fas']=length(intersect(gtex$id,sub))# takes a while
}
sens=sens/length(gtex$id)

pdf(sprintf('%s/eAssociations.ROC.pdf',fig_dir))
plot(FDR,sens[,'ras'],type='l',col='black',xlab='FDR',ylab='sensitivity')
lines(FDR,sens[,'fas'],col='red')
legend('topleft',col=c('black','red'),lty=1,legend=c('ras','fas'))
dev.off()


# Overlap FastQTL and RASQUAL:
merged=merge(fas[,list(SNP,gene,p.value)],ras[,list(SNP,gene,p.value)],,by=c('SNP','gene'),all=TRUE,suffixes = c(".fas",".ras"))
c=merged[,list(ras=!is.na(p.value.ras),fas=!is.na(p.value.fas))]

pdf(sprintf('%s/eAssociations.venn.pdf',fig_dir))
a=vennCounts(c)
vennDiagram(a)
a[,'Counts']=round(a[,'Counts']/sum(a[,'Counts']),4)
vennDiagram(a)
dev.off()