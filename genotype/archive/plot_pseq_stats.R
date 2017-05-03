#!/usr/bin/env Rscript
source('utils.R')

plot_istats = function(istats){
	# fields description: 
	# ID 		Individual ID
	# NALT		Number of non-reference genotypes
	# NMIN		Number of genotypes with a minor allele
	# NHET		Number of heterozygous genotypes for individual
	# NVAR		Total number of called variants for individual
	# RATE		Genotyping rate for individual
	# SING		Number of singletons individual has
	# TITV		Mean Ti/Tv for variants for which individual has a nonreference genotype
	# PASS		Number of variants PASS'ing for which individual has a nonreference genotype
	# PASS_S	Number of singletons PASS'ing for which individual has a (singleton) nonreference genotype
	# QUAL		Mean QUAL for variants for which individual has a nonreference genotype
	# DP		Mean variant DP for variants for which individual has a nonreference genotype
	x = names(istats)[1]
	ys = names(istats)[2:ncol(istats)]
	names(ys) = c('# non-reference genotypes','# genotypes with a minor allele', '# heterozygous genotypes', '# called variants', 'Genotyping rate', 'Number of singletons','Mean Ti/Tv for variants with nonreference genotype','# PASSes for variants with nonreference genotype', '# PASSes for variants with nonreference genotype', 'Mean QUAL for variants with nonreference genotype', 'Mean depth for variants with nonreference genotype')
	p = list()
	for (i in seq(1,10)){
		p[[i]] =  ggplot(istats, aes_string(x, ys[i], label = ys[i], group = 1)) + geom_line() + xlab(x) + ylab(names(ys)[i]) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + background_grid(major = "xy", minor = "xy") + geom_text(angle = 90, vjust = 0.5, hjust = 0)
	}
	p_grid = plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]], ncol = 2)
	return(p_grid)
}

#----- main -------
args = commandArgs(TRUE)
cw = args[1]
setwd(cw)

# before genotype refinement: 
istats_path = 'recalibrated_variants.vcf.pseq.istats'
istats1 = fread(istats_path)
istats1$procedure = 'before refinement'
figure_path = '../figures/pseq_istats_recalibrated_variants.pdf'
p_grid = plot_istats(istats1)
save_plot(figure_path, p_grid, base_height = 40, base_width = 32)



# after genotype refinement: 
istats_path = 'recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf.pseq.istats'
istats2 = fread(istats_path)
istats2$procedure = 'after refinement'
figure_path = '../figures/pseq_istats_recalibrated_variants.GRCh37.postCGP.GQfiltered.pdf'
p_grid = plot_istats(istats2)
save_plot(figure_path, p_grid, base_height = 40, base_width = 32)



# before and after overlay: 
istats_rbind = rbind(istats1, istats2) 
x = names(istats_rbind)[1]
ys = names(istats_rbind)[2:(ncol(istats_rbind)-1)]
names(ys) = c('# non-reference genotypes','# genotypes with a minor allele', '# heterozygous genotypes', '# called variants', 'Genotyping rate', 'Number of singletons','Mean Ti/Tv for variants with nonreference genotype','# PASSes for variants with nonreference genotype', '# PASSes for variants with nonreference genotype', 'Mean QUAL for variants with nonreference genotype', 'Mean depth for variants with nonreference genotype')
p = list()
for (i in seq(1,10)){
	p[[i]] =  ggplot(istats_rbind, aes_string(x, ys[i], group = 'procedure', color = 'procedure')) + geom_line() + xlab(x) + ylab(names(ys)[i]) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + background_grid(major = "xy", minor = "xy")
}
p_grid = plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]], ncol = 2)
figure_path = '../figures/pseq_istats_overlay.pdf'
save_plot(figure_path, p_grid, base_height = 40, base_width = 32)
