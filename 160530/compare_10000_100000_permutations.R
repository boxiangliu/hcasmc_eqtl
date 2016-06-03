#!/usr/bin/env Rscript
# bosh liu
# durga
# compare 10000 and 100000 permutations

# command args:
args=commandArgs(T)
perm10000_file=args[1]
perm100000_file=args[2]
figure=args[3]

# read input: 
perm10000=fread(perm10000_file,header=F)
perm100000=fread(perm100000_file,header=F)


# setnames:
setnames(perm10000, c('gene_id','num_variants','shape1','shape2','dummy','best_variant','distance','nominal_p','slope','permute_p','beta_p'))
setnames(perm100000, c('gene_id','num_variants','shape1','shape2','dummy','best_variant','distance','nominal_p','slope','permute_p','beta_p'))


# make scatterplot: 
pdf(figure)
plot(perm10000$beta_p,perm100000$beta_p, xlab="max 10000 permutations",ylab='max 100000 permutations')
abline(0,1)
dev.off()

