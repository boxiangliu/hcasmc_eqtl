# plot the number of significant sqtl vs distance to intron start site. 
library(cowplot)
library(stringr)
library(data.table)
library(dplyr)
library(dtplyr)

# command args:
args=commandArgs(T)
sqtl_file=args[1]
fig_dir=args[2]
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}


# read sqtl:
if (str_detect(sqtl_file,'.gz')){
	sqtl=fread(sprintf('zcat %s',sqtl_file),col.names=c('pheno','snp','dist','pval','beta','se'))
} else {
	sqtl=fread(sqtl_file,col.names=c('pheno','snp','dist','pval','beta','se'))
}


# subset to significant sqtls:
sig_sqtl=sqtl[pval<1e-4,]


# plot number of sig. sqtls vs distance:
p=ggplot(sig_sqtl,aes(dist))+geom_histogram(binwidth=100)
save_plot(sprintf('%s/sqtl_vs_dist.pdf',fig_dir),p)

# calculate length of intron:
temp=str_split_fixed(sig_sqtl$pheno,":",4)
sig_sqtl$start=as.numeric(temp[,2])
sig_sqtl$end=as.numeric(temp[,3])
sig_sqtl=sig_sqtl%>%mutate(length=end-start)


# calculate distance to intron start as a fraction of intron length:
sig_sqtl=sig_sqtl%>%mutate(fraction=dist/length)
sig_sqtl_within_intron=sig_sqtl[fraction<=1&fraction>=0]


# plot number of sig. sqtls vs fractional distance from intron start: 
p=ggplot(sig_sqtl_within_intron,aes(fraction))+geom_histogram(binwidth=0.01)+xlab('fraction within intron')
save_plot(sprintf('%s/intron_sqtl.pdf',fig_dir),p)
