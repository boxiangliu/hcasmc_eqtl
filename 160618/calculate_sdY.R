#!/usr/bin/env Rscript
# bosh liu
# durga
# plot the expression vs genotype


# command args:
args=commandArgs(T,T)
expression_file=args$expression
covariate_file=args$covariate
output_file=args$output
# expression_file='../processed_data/160530/combined.filter.norm.rpkm'
# covariate_file='../processed_data/160530/find_optimum_num_PEER_factors_matrixeqtl/covariates.matrixeqtl.pc3.peer8.tsv'


# read input:
expression=fread(expression_file,header=T)
covariate=read.table(covariate_file,header=T,row.names=1,check.name=F,stringsAsFactors=F)


# cast data.table into data.frame, and
# make the first column the row names:
expression=as.data.frame(expression)
rownames(expression)=expression[,1]
expression=expression[,-1]


# cast covariate into a design matrix:
covariate['gender',]=ifelse(covariate['gender',]=='M',0,1)
covariate=as.matrix(t(covariate))
covariate[,'gender']=covariate[,'gender']-1


sd_df=data.frame()
counter=0
for (gene_id in rownames(expression)){
	counter=counter+1
	if (counter%%1000==0){
		message(counter,' genes processed.')
	}


	# regress covariate out from expression:
	gene=t(expression[rownames(expression)==gene_id,])
	stopifnot(names(gene)==rownames(covariate))
	fit=lm(gene~covariate)
	rsd=fit$residuals


	# calculate standard deviation:
	sd=sd(rsd)
	sd_df=rbind(sd_df,data.frame(gene_id=gene_id,sd=sd))
}


# write output:
write.table(sd_df,output_file,sep='\t',col.names=T,row.names=F,quote=F)
