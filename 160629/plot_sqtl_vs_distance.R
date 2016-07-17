# plot the number of significant sqtl vs distance to intron start site. 


# command args:
args=commandArgs(T,T)
sqtl_file=args$input
all_sqtl_fig=args$all_sqtl_fig
intronic_sqtl_fig=args$intronic_sqtl_fig

# sqtl_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160629//sqtl.nominal.allpairs.normal.1e5.padj.txt'
# all_sqtl_fig='/srv/persistent/bliu2/HCASMC_eQTL/figures/160629/num_sig_sqtl_vs_dist.pdf'
# intronic_sqtl_fig='/srv/persistent/bliu2/HCASMC_eQTL/figures/160629/num_sig_sqtl_within_intron_vs_dist.pdf'

# read sqtl:
sqtl=fread(sqtl_file)


# subset to significant sqtls:
sig_sqtl=sqtl[qval<0.05,]


# plot number of sig. sqtls vs distance:
p=ggplot(sig_sqtl,aes(dist))+geom_histogram(binwidth=100)
save_plot(all_sqtl_fig,p)

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
save_plot(intronic_sqtl_fig,p)
