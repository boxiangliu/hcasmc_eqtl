library(data.table)
library(stringr)
library(foreach)
library(doMC)
registerDoMC(10)
library(cowplot)

in_dir = '../data/atacseq/fbs/'
fig_dir = '../figures/atacseq/qc/plot_peak_count/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

read_narrowPeak = function(fn){
	command = sprintf('zcat %s',fn)
	narrowPeak = fread(command)
	setnames(narrowPeak,c('chrom','chromStart','chromEnd','name','score','strand','signalValue','pValue','qValue','peak'))
	return(narrowPeak)
}

count_peak_by_signalValue = function(narrowPeak, SIGNALVALUE){
	sum(narrowPeak$signalValue > SIGNALVALUE)
}

count_peak_by_FDR = function(narrowPeak,FDR){
	FDR = -log10(FDR)
	sum(narrowPeak$qValue > FDR)
}

plot_peak_count = function(count){
	melted_count = melt(count,id.vars = 'sample')
	setorder(melted_count,-value)
	melted_count[,sample := factor(sample,unique(sample))]
	melted_count[,variable := factor(variable,c('n','n_by_FDR','n_by_signalValue'))]
	ggplot(melted_count, aes(sample,value,fill=variable)) + 
		geom_bar(stat = 'identity', position = 'dodge') + 
		scale_y_log10() + 
		xlab('') + 
		ylab('Number of peaks') + 
		scale_fill_discrete(name = '', breaks = c('n_by_FDR','n_by_signalValue'), labels = c('FDR < 0.05', 'Fold Enrichment > 5')) +
		coord_flip() 
}

fn_list = list.files(in_dir,pattern='sorted.narrowPeak.gz',full.names=TRUE,recursive=TRUE)
fn_list = fn_list[fn_list != "../data/atacseq/fbs//2108/out/peak/macs2/overlap/2108.sorted.narrowPeak.gz"]
count = foreach(fn = fn_list, .combine = rbind)%dopar%{
	narrowPeak = read_narrowPeak(fn)
	sample = basename(fn)
	sample = str_replace(sample,'.sorted.narrowPeak.gz','')
	n_by_signalValue = count_peak_by_signalValue(narrowPeak,5)
	n_by_FDR = count_peak_by_FDR(narrowPeak,0.05)
	data.table(n_by_signalValue,n_by_FDR,sample)
}

p = plot_peak_count(count) + theme(legend.position = 'top')
fig_fn = sprintf('%s/peak_count.pdf',fig_dir)
save_plot(fig_fn, p, base_width = 4, base_height = 4)

summary(count$n_by_FDR)
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 # 109533  177415  289614  250521  309210  350303 