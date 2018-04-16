library(data.table)
library(stringr)
library(foreach)
library(cowplot)

samtools_stats_dir = '../processed_data/atacseq/qc/samtools_stats/'
fig_dir = '../figures/atacseq/qc/plot_samtools_stats/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

read_samtools_stats = function(fn){
	command = sprintf('cat %s | grep ^SN | cut -f 2- ', fn)
	sample = basename(fn)
	stats = fread(command,sep=':')
	setnames(stats,c('field','value'))
	stats[,value := trimws(value)]
	stats[,value := str_split_fixed(value,'\\t',2)[,1]]
	stats[,value := as.numeric(value)]
	stats$sample = sample
	return(stats)
}

munge_sample_name = function(stats){
	stats[sample == 'CA1508_L1_TAAGGCGA_L001_R1.trim.PE2SE.bam.stats', sample := '1508-rep1']
	stats[sample == 'CA1508_L2_TAAGGCGA_L002_R1.trim.PE2SE.bam.stats', sample := '1508-rep2']
	stats[,sample:=str_split_fixed(sample,'_',2)[,1]]
	stats[,sample:=str_replace(sample,'CA','')]
	stats[,sample:=str_replace(sample,'FBS','rep')]
	return(stats)
}

plot_read_depth = function(stats){
	x = stats[field %in% c('raw total sequences','reads mapped')]
	setorder(x,-value)
	x[,sample:=factor(sample,unique(sample))]
	ggplot(x, aes(sample,round(value/2),fill=field)) + 
		geom_bar(stat = 'identity', position = 'dodge') +
		scale_y_log10() + 
		xlab('') +
		ylab('Number of reads') + 
		scale_fill_discrete(name = '', breaks = c('raw total sequences','reads mapped'), labels = c('Total reads', 'Mapped reads')) +
		coord_flip()
}

fn = '../processed_data/atacseq/qc/samtools_stats/CA1508_L2_TAAGGCGA_L002_R1.trim.PE2SE.bam.stats'
fn_list = list.files(samtools_stats_dir,pattern = 'stats', full.names=TRUE)
stats = foreach(fn = fn_list,.combine = rbind)%do%{
	read_samtools_stats(fn)
}


stats = munge_sample_name(stats)
stats = stats[sample!='2108']
p = plot_read_depth(stats) + theme(legend.position="top")
fig_fn = sprintf('%s/number_of_reads.pdf',fig_dir)
save_plot(fig_fn,p,base_width = 4, base_height = 4)
summary(stats[field == 'raw total sequences',round(value/2)])
 #     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
 # 18042227  37705038  44604463  57728893  69777832 137248057summary(stats[field == 'reads mapped',round(value/2)])
stats[field == 'maximum length']
#              field value    sample
#  1: maximum length    50    200212
#  2: maximum length    50      1346
#  3: maximum length    50 1508-rep1
#  4: maximum length    50 1508-rep2
#  5: maximum length    76      1522
#  6: maximum length    76      2108
#  7: maximum length    76 2305-rep1
#  8: maximum length    76 2305-rep2
#  9: maximum length    76 2305-rep3
# 10: maximum length    76      2356
# 11: maximum length    76      2510
# 12: maximum length    76      2989
summary(stats[field == 'reads mapped',round(value/2)])
 #     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
 # 17786536  33349461  37166884  52967857  67776813 135123290