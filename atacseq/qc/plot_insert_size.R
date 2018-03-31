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

data = foreach(i = seq_along(file_list),.combine='rbind')%do%{
	sample_name = names(file_list)[i]
	dt = read_insert_size_log_file(file_list[i])
	setnames(dt,c('insert_size','count'))
	dt$sample_name = sample_name
	return(dt)
}

p = plot_insert_size_distr(data)
fig_path = sprintf('%s/insert_size_distr.pdf',fig_dir)
save_plot(fig_path, p, base_width = 6, base_height = 4)