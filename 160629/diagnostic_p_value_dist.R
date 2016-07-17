# diagnose the enrichment of high p-values:

# read leafcutter counts
count=fread('zcat /srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/leafcutter/leafcutter_perind_numers.counts.gz',header=F)%>%as.data.frame()
rownames(count)=count[,1]
count=count[,-1]
count=as.matrix(count)


# calculate mean count, the number of samples with zero reads, the number of samples with greater than 6 reads. 
mean_count=rowMeans(count)
mean_count=data.frame(pheno=names(mean_count),mean_count)
zero_count=apply(count,1,function(x) sum(x==0))
zero_count=data.frame(pheno=names(zero_count),zero_count)
gt6_count=apply(count,1,function(x) sum(x>6))
gt6_count=data.frame(pheno=names(gt6_count),gt6_count)


# read sqtl file:
sqtl=fread('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160629//sqtl.nominal.allpairs.normal.1e5.txt')
setnames(sqtl, c('pheno','geno','dist','pval','beta','varbeta'))


# calculate mean p-value, minimum p-value and adjusted minimum p-value:
# adjusted p-value = p-value*number of tests
mean_pval=sqtl%>%group_by(pheno)%>%summarize(mean_pval=mean(pval))
min_pval=sqtl%>%group_by(pheno)%>%summarize(min_pval=min(pval),n=n())
min_pval=min_pval%>%mutate(adj_min_pval=ifelse(min_pval*n>1,1,min_pval*n))


# merge all relavant data frames:
merged=merge(mean_count,mean_pval,by='pheno')
merged=merge(merged,min_pval,by='pheno')
merged=merge(merged,zero_count,by='pheno')
merged=merge(merged,gt6_count,by='pheno')


# plotting
# mean count vs mean pvalue:
p=ggplot(merged,aes(mean_count,mean_pval))+stat_binhex()+scale_x_log10()
save_plot('/srv/persistent/bliu2/HCASMC_eQTL/figures/160629/mean_count_vs_mean_pval.pdf',p)

# mean count vs mean pvalue boxplot:
merged=merged%>%mutate(log_mean_count=log(mean_count+2))
merged=merged%>%mutate(log_mean_count_bin=cut(log_mean_count,breaks=30,include.lowest=T))
p=ggplot(merged,aes(log_mean_count_bin,mean_pval))+geom_boxplot()+geom_hline(yintercept=0.5)+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
save_plot('/srv/persistent/bliu2/HCASMC_eQTL/figures/160629/mean_count_vs_mean_pval.boxplot.pdf',p)


# mean count vs minimum pvalue: 
p=ggplot(merged,aes(mean_count,min_pval))+stat_binhex()+scale_x_log10()
save_plot('/srv/persistent/bliu2/HCASMC_eQTL/figures/160629/mean_count_vs_min_pval.pdf',p)


# mean count vs adjusted minimum pvalue: 
p=ggplot(merged,aes(mean_count,adj_min_pval))+stat_binhex()+scale_x_log10()
save_plot('/srv/persistent/bliu2/HCASMC_eQTL/figures/160629/mean_count_vs_adj_min_pval.pdf',p)


# zero count vs mean pvalue: 
p=ggplot(merged,aes(zero_count,mean_pval))+stat_binhex()
save_plot('/srv/persistent/bliu2/HCASMC_eQTL/figures/160629/zero_count_vs_mean_pval.pdf',p)


# zero count vs minimum pvalue:
p=ggplot(merged,aes(zero_count,min_pval))+stat_binhex()
save_plot('/srv/persistent/bliu2/HCASMC_eQTL/figures/160629/zero_count_vs_min_pval.pdf',p)


# zero count vs adjusted mean pvalue: 
p=ggplot(merged,aes(zero_count,adj_min_pval))+stat_binhex()
save_plot('/srv/persistent/bliu2/HCASMC_eQTL/figures/160629/zero_count_vs_adj_mean_pval.pdf',p)


# greater than 6 count vs mean pvalue: 
p=ggplot(merged,aes(gt6_count,mean_pval))+stat_binhex()
save_plot('/srv/persistent/bliu2/HCASMC_eQTL/figures/160629/gt6_count_vs_mean_pval.pdf',p)


# greater than 6 count vs mean pvalue boxplot: 
p=ggplot(merged,aes(as.factor(gt6_count),mean_pval))+geom_boxplot()+geom_hline(yintercept=0.5,color='red')
save_plot('/srv/persistent/bliu2/HCASMC_eQTL/figures/160629/gt6_count_vs_mean_pval.boxplot.pdf',p)


# greater than 6 count vs mean pvalue: 
p=ggplot(merged,aes(gt6_count,min_pval))+stat_binhex()
save_plot('/srv/persistent/bliu2/HCASMC_eQTL/figures/160629/gt6_count_vs_min_pval.pdf',p)


# greater than 6 count vs mean pvalue: 
p=ggplot(merged,aes(gt6_count,adj_min_pval))+stat_binhex()
save_plot('/srv/persistent/bliu2/HCASMC_eQTL/figures/160629/gt6_count_vs_adj_min_pval.pdf',p)
