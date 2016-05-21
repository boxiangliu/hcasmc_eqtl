#!/usr/bin/env Rscript
# bosh liu
# durga
# 16/05/16
# to study whether omitting covariates will reduce power

# library:
library(gap)


# paths:
figure_path='../figures/160516_sim_study_on_covariates/'


# constants:
num_genes=20
num_snps=20
num_samples=50


# simulate sex:
sex=sample(c(0,1),size=num_samples,replace=T)


# simulate age:
age=sample(seq(10,90),size=num_samples,replace=T)


# simulate ethinicity: 
race=sample(seq(1,4),size=num_samples,replace=T)


# simulate genotypes:
genotypes=sample(seq(0,2),size=num_snps*num_samples,replace=T)


# simulate gene expression:
genotype_eff=5
expression=rnorm(num_snps*num_samples,mean=genotype_eff*genotypes,sd=1)


# cast genotype and expression into matrices:
genotypes_mat=matrix(genotypes,nrow=num_snps,ncol=num_samples)
expression_mat=matrix(expression,nrow=num_genes,ncol=num_samples)


# add covariates:
sex_eff=10
sex_cov=rnorm(length(sex),sex_eff*sex,sd=1)
age_eff=0.1
age_cov=rnorm(length(age),age_eff*age,sd=1)
race_eff=2.5
race_cov=rnorm(length(race),race_eff*race,sd=1)
expression_cov=t(t(expression_mat)+sex_cov+age_cov+race_cov+6)


# eQTL mapping:
pval_list1=pval_list2=pval_list3=pval_list4=c()
for (i in seq(1,num_snps)){
	for (j in seq(1,num_genes)){
		res=summary(lm(expression_cov[i,]~genotypes_mat[j,]+sex_cov+age_cov+race_cov))
		pval_list1=c(pval_list1,coef(res)[2,4])
		
		res=summary(lm(expression_cov[i,]~genotypes_mat[j,]+sex_cov+age_cov))
		pval_list2=c(pval_list2,coef(res)[2,4])
		
		res=summary(lm(expression_cov[i,]~genotypes_mat[j,]+sex_cov))
		pval_list3=c(pval_list3,coef(res)[2,4])
		
		res=summary(lm(expression_cov[i,]~genotypes_mat[j,]))
		pval_list4=c(pval_list4,coef(res)[2,4])
	}
}


# qqplot:
pdf(paste0(figure_path,'qqplots.pdf'))
qqunif(pval_list1)
points(qqunif(pval_list2,plot=F),col='red')
points(qqunif(pval_list3,plot=F),col='green')
points(qqunif(pval_list4,plot=F),col='black')
legend('topleft',legend=c('~genotype+sex+age+race','~genotype+sex+age','~genotype+sex','~genotype'),col=c('blue','red','green','black'),pch=1)
dev.off()