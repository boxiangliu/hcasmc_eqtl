#!/usr/bin/env Rscript
# bosh liu
# durga
# plot the expression vs genotype

# library:
library(cowplot)
library(data.table)


# command args:
args=commandArgs(T)
expression_file=args[1]
genotype_file=args[2]
covariate_file=args[3]
# gene_id=args[4]
# snp_id=args[5]
fastqtl_file=args[4]
figure_dir=args[5]
# expression_file='../processed_data/160530/combined.filter.norm.rpkm'
# genotype_file='../processed_data/160530/dosage.tsv'
# covariate_file='../processed_data/160530/covariates.tsv'
# fastqtl_file='../processed_data/160530/fastqtl.padj.txt'



# read input:
expression=fread(expression_file,header=T)
genotype=fread(genotype_file,header=T)
covariate=read.table(covariate_file,header=T,row.names=1,check.name=F,stringsAsFactors=F)
fastqtl=fread(fastqtl_file,header=T)


# cast data.table into data.frame, and
# make the first column the row names:
expression=as.data.frame(expression)
genotype=as.data.frame(genotype)
rownames(expression)=expression[,1]
rownames(genotype)=genotype[,1]
expression=expression[,-1]
genotype=genotype[,-1]


# cast covariate into a design matrix:
# covariate['gender',]=ifelse(covariate['gender',]=='M',0,1)
covariate=as.matrix(t(covariate))
covariate[,'gender']=covariate[,'gender']-1


# sort the fastqtl file by p-value:
rk=rank(fastqtl$beta_p,ties.method='min')
ord=order(fastqtl$beta_p)
fastqtl$rank=rk
fastqtl=fastqtl[ord,]


for (i in seq(nrow(fastqtl))){
	gene_id=fastqtl[i,gene_id]
	snp_id=fastqtl[i,best_variant]
	p_value=fastqtl[i,beta_p]
	q_value=fastqtl[i,qvalues]
	rk=fastqtl[i,rank]
	figure_path=paste0(figure_dir,"/",sprintf("%05d",rk),"_",gene_id,"_",snp_id,".pdf")
	message('plotting ',figure_path)

	# regress covariate out from expression:
	gene=t(expression[rownames(expression)==gene_id,])
	stopifnot(names(gene)==rownames(covariate))
	fit=lm(gene~covariate)
	rsd=fit$residuals


	# get genotype: 
	snp=unlist(genotype[rownames(genotype)==snp_id,])
	snp=round(snp)


	# plot genotype vs expression: 
	to_plot=data.frame(snp,rsd,gene)

	p=ggplot(to_plot,aes(x=as.factor(snp),y=rsd))+
		geom_boxplot(outlier.size=-1)+
		geom_jitter(width=0.1)+
		stat_summary(fun.y=mean, colour="darkred", geom="point", 
	               shape=18, size=5,show.legend = FALSE) + 
		ggtitle(paste0('p-value = ',signif(p_value,4),'\n','q-value = ',signif(q_value,4))) +
		xlab(snp_id)+ylab(gene_id)

	save_plot(figure_path,p)
}

