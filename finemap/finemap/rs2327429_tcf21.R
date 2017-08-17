library(data.table)
library(stringr)
library(cowplot)

fig_dir='../figures/finemap/finemap/rs2327429_tcf21/'
if (!dir.exists(fig_dir)) dir.create(fig_dir,recursive=TRUE)

# Genotype: 
cmd="bcftools query -H -r chr6:134209837-134209837 -f '%CHROM\t%POS\t%ID[\t%GT]\n' ../data/joint3/asvcf/phased_and_imputed.chr6.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz > temp"
system(cmd)
rs2327429=fread('temp')
unlink('temp')

colnames(rs2327429)=
	str_split_fixed(colnames(rs2327429),']|:',3)[,2]

rs2327429=melt(rs2327429,
	id.vars=c('ID','CHROM','POS'),
	variable.name='sample',
	value.name='genotype')

hap_to_geno=c('0|0'=0,'0|1'=1,'1|0'=1,'1|1'=2)
rs2327429[,genotype:=hap_to_geno[genotype]]


# Expression:
expression=fread('../data/rnaseq2/read_count/rnaseqc/rnaseqc.hcasmc_eqtl.reads.gct',header=TRUE)
expression=expression[Description%in%c('TCF21','RP3-323P13.2'),]
temp=melt(expression,id.vars=c('Name','Description'),variable.name='sample')
expression=dcast(temp,sample~Description)


# Covariates: 
cmd='Rscript rasqual/combine_covariates.R --genotype_pc=../processed_data/genotype/genotype_pc/genotype_pcs.52samples.tsv --peer=../processed_data/eqtl/peer/factors.tsv --sample_info=../data/sample_info/sample_info.xlsx --output=temp --gender_coding=numerical --num_geno_pc=4 --num_peer_factor=8'
system(cmd)
covariate=read.table('temp',col.names=c(paste0('PC',1:4),paste0('PEER',1:8),'sex'))
rownames(covariate)=expression$sample
covariate$sample=rownames(covariate)
setDT(covariate)
unlink('temp')


# Offset:
offset=fread('../processed_data/rasqual/expression/K.txt',col.names=c('Name',covariate$sample))
offset=offset[Name=='ENSG00000118526.6']
offset=melt(offset,id.var='Name',variable.name='sample',value.name='offset')

pdf(sprintf('%s/tcf21_vs_offset.pdf',fig_dir))
plot(offset$offset,expression$TCF21,log='y',ylab='TCF21',xlab='Offset')
dev.off()


# Take residual:
expression[,c('TCF21','RP3-323P13.2'):=list(log2(TCF21+1),log2(`RP3-323P13.2`+1))]
merged=merge(expression,offset,by='sample')
merged=merge(merged,covariate,by='sample')
fit=lm(TCF21~PC1+PC2+PC3+PC4+PEER1+PEER2+PEER3+PEER4+PEER5+PEER6+PEER7+PEER8+sex+offset,data=merged)
summary(fit)
expression$TCF21_residuals=fit$residuals

fit=lm(`RP3-323P13.2`~PC1+PC2+PC3+PC4+PEER1+PEER2+PEER3+PEER4+PEER5+PEER6+PEER7+PEER8+sex+offset,data=merged)
summary(fit)
expression$TARID_residuals=fit$residuals

# Plot:
merged=merge(rs2327429,expression[,list(sample,TCF21_residuals,TARID_residuals)],by='sample')
p1=ggplot(merged,aes(as.character(genotype),TCF21_residuals))+geom_boxplot(outlier.size=-1)+geom_jitter(width=0.1)+xlab('Genotype')+ylab('Expression residual')+scale_x_discrete(labels=c('0'='TT','1'='TC','2'='CC'))+ggtitle('TCF21:rs2327429')
p2=ggplot(merged,aes(as.character(genotype),TARID_residuals))+geom_boxplot(outlier.size=-1)+geom_jitter(width=0.1)+xlab('Genotype')+ylab('Expression residual')+scale_x_discrete(labels=c('0'='TT','1'='TC','2'='CC'))+ggtitle('TARID:rs2327429')
pdf(sprintf('%s/residual_vs_genotype.pdf',fig_dir))
p1;p2
dev.off()