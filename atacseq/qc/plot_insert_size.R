library(stringr)
library(data.table)
library(foreach)
library(cowplot)

# Variables:
atacseq_dir = '../data/atacseq/fbs/'
fn_pattern = 'inserts.hist_data.log'
fig_dir = '../figures/atacseq/qc/plot_insert_size/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

# Functions:
get_insert_size_log_file = function(dir,pattern){
	file_list = list.files(dir,pattern,full.names=TRUE,recursive=TRUE)
	replicate = basename(dirname(file_list))
	sample_name = basename(file_list)
	sample_name = str_extract(sample_name,'(?<=^)(.+?)(?=[-_])') 
	sample_name = str_replace(sample_name,'CA','')
	sample_name = paste(sample_name,replicate,sep=' ')	
	names(file_list) = sample_name
	return(file_list)
}


read_insert_size_log_file = function(fn){
	dt = fread(fn,skip=10)
	return(dt)
}


plot_insert_size_distr = function(data){
	p = ggplot(data,aes(x = insert_size, y = count, fill = sample_name))+
		geom_area()+
		scale_fill_discrete(name = 'Sample')+
		xlab('Insert size') + 
		ylab('Read count')
	return(p)
}

# Main: 
file_list = get_insert_size_log_file(atacseq_dir,fn_pattern)
file_list = file_list[file_list != "../data/atacseq/fbs//2108/out/qc/rep1/CA2108_S4_concat_R1_001.PE2SE.inserts.hist_data.log"]

data = foreach(i = seq_along(file_list),.combine='rbind')%do%{
	sample_name = names(file_list)[i]
	dt = read_insert_size_log_file(file_list[i])
	setnames(dt,c('insert_size','count'))
	dt$sample_name = sample_name
	return(dt)
}

p = plot_insert_size_distr(data) + 
	scale_fill_discrete(name = 'Sample', breaks = sort(c(paste0('1508 rep',1:2),paste0('2305 rep',1:3),paste0(c('1346','1522','2356','2510','2989'),' rep1'))), labels = sort(c(paste0('1508-rep',1:2),paste0('2305-rep',1:3),c('1346','1522','2356','2510','2989')))) +
	theme(legend.position = c(0.95, 0.95),legend.justification = c("right", "top"))
fig_path = sprintf('%s/insert_size_distr.pdf',fig_dir)
save_plot(fig_path, p, base_width = 4, base_height = 4)