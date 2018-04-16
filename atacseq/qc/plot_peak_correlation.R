library(data.table)
library(stringr)
library(foreach)
library(doMC)
registerDoMC(10)
library(cowplot)
library(ggcorrplot)
in_dir = '../data/atacseq/fbs/'
fig_dir = '../figures/atacseq/qc/plot_peak_correlation/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

calculate_jaccard = function(a,b,sorted = TRUE){
	if (sorted == FALSE){
		sorted_a = paste0(a,'.sorted.bed')
		on.exit(unlink(sorted_a))
		sorted_b = paste0(b,'.sorted.bed')
		on.exit(unlink(sorted_b))
		command = sprintf('zcat %s | sort -k1,1 -k2,2n > %s',a,sorted_a)
		system(command)
		command = sprintf('zcat %s | sort -k1,1 -k2,2n > %s',b,sorted_b)
		system(command)
		command = sprintf('bedtools jaccard -a %s -b %s',sorted_a,sorted_b)
		result = system(command,intern=TRUE)
	} else {
		command = sprintf('bedtools jaccard -a %s -b %s',a,b)
		result = system(command,intern=TRUE)
	}
	return(result)
}

make_correlation_matrix = function(jaccard){
	jaccard_flip = copy(jaccard)
	setnames(jaccard_flip,c('sample1','sample2'),c('sample2','sample1'))
	all_sample = unique(c(jaccard$sample1,jaccard$sample2))
	self = data.table(sample1 = all_sample, sample2 = all_sample, jaccard = 1)
	jaccard = rbind(jaccard,jaccard_flip,self)
	jaccard = dcast(jaccard,sample1~sample2,value.var='jaccard')
	setDF(jaccard)
	rownames(jaccard) = jaccard[,'sample1']
	jaccard[,'sample1'] = NULL
	return(jaccard)
}

#---------------#
# Cross samples #
#---------------#
fn_list = list.files(in_dir,pattern='sorted.narrowPeak.gz',full.names=TRUE,recursive=TRUE)
pairs = combn(fn_list,2)

jaccard = foreach(i = 1:ncol(pairs), .combine = rbind)%dopar%{
	a = pairs[1,i]
	b = pairs[2,i]
	sample1 = basename(a)
	sample1 = str_replace(sample1,'.sorted.narrowPeak.gz','') 
	sample2 = basename(b)
	sample2 = str_replace(sample2,'.sorted.narrowPeak.gz','') 
	jaccard = calculate_jaccard(a,b)
	colnames = jaccard[1]
	jaccard = jaccard[2]
	colnames = str_split_fixed(colnames,'\\t',4)
	jaccard = str_split_fixed(jaccard,'\\t',4)
	jaccard = data.table(jaccard[1],jaccard[2],jaccard[3],jaccard[4])
	setnames(jaccard,colnames)
	jaccard$sample1 = sample1
	jaccard$sample2 = sample2
	return(jaccard)
}

jaccard = jaccard[, list(jaccard = as.numeric(jaccard),sample1,sample2)]
jaccard_cor = make_correlation_matrix(jaccard)
colnames(jaccard_cor) = paste0('S',colnames(jaccard_cor))
rownames(jaccard_cor) = paste0('S',rownames(jaccard_cor))
p = ggcorrplot(as.matrix(jaccard_cor),type = 'lower',lab = TRUE)
fig_fn = sprintf('%s/cross_sample_correlation.pdf',fig_dir)
save_plot(fig_fn,p)

#------# 
# 2305 #
#------#
fn_list = c("../data/atacseq/fbs//2305/out/peak/macs2//rep1/CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.nodup.tn5.pf.pval0.1.narrowPeak.gz",
		"../data/atacseq/fbs//2305/out/peak/macs2//rep2/CA2305-FBS2_S2_concat_R1_001.trim.PE2SE.nodup.tn5.pf.pval0.1.narrowPeak.gz",
		"../data/atacseq/fbs//2305/out/peak/macs2//rep3/CA2305-FBS3_S3_concat_R1_001.trim.PE2SE.nodup.tn5.pf.pval0.1.narrowPeak.gz")
pairs = combn(fn_list,2)
jaccard = foreach(i = 1:ncol(pairs), .combine = rbind)%dopar%{
	a = pairs[1,i]
	b = pairs[2,i]
	sample1 = basename(a)
	# sample1 = str_replace(sample1,'.sorted.narrowPeak.gz','') 
	sample1 = str_split_fixed(sample1,'_',2)[,1]
	sample2 = basename(b)
	# sample2 = str_replace(sample2,'.sorted.narrowPeak.gz','') 
	sample2 = str_split_fixed(sample2,'_',2)[,1]
	jaccard = calculate_jaccard(a,b,sorted=FALSE)
	colnames = jaccard[1]
	jaccard = jaccard[2]
	colnames = str_split_fixed(colnames,'\\t',4)
	jaccard = str_split_fixed(jaccard,'\\t',4)
	jaccard = data.table(jaccard[1],jaccard[2],jaccard[3],jaccard[4])
	setnames(jaccard,colnames)
	jaccard$sample1 = sample1
	jaccard$sample2 = sample2
	return(jaccard)
}

jaccard = jaccard[, list(jaccard = as.numeric(jaccard),sample1,sample2)]
jaccard_cor = make_correlation_matrix(jaccard)
p = ggcorrplot(as.matrix(jaccard_cor))
fig_fn = sprintf('%s/2305_rep_correlation.pdf',fig_dir)
save_plot(fig_fn,p)