library(stringr)
library(data.table)
library(foreach)
library(cowplot)
library(heatmap3)

in_dir = '../processed_data/atacseq/qc/plot_vplot/'
pattern = 'profile.txt'
matrix_fn = '../processed_data/atacseq/qc/plot_vplot/1346_rep1.matrix.gz'
fig_dir = '../figures/atacseq/qc/plot_vplot.R/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir)}

get_profile_file_list = function(dir,pattern){
	file_list = list.files(in_dir,pattern=pattern,full.names=TRUE,recursive=TRUE)
	sample_name = basename(file_list)
	sample_name = str_extract(sample_name,'(?<=^)(.+?)(?=\\.profile)')
	sample_name = str_replace(sample_name,'_',' ')
	names(file_list) = sample_name
	return(file_list)
}

read_matrix = function(fn){
	df = read.table(fn,comment.char='@')
	row_data = df[,1:6]
	names(row_data) = c('chr','start','end','gene','unknown','strand')
	df = df[,7:ncol(df)]
	return(list(row_data,df))
}

read_profile = function(fn,sample_name){
	dt = fread(fn)
	dt = dt[,3:ncol(dt)]
	dt = data.table(t(dt))
	setnames(dt,c('bin','insertions'))
	dt$sample_name = sample_name
	return(dt)
}

plot_profile = function(data){
	p = ggplot(data,aes(x = bin, y = insertions, color = sample_name)) + 
		geom_line() + 
		xlab('Bin') + 
		ylab('Insertion') + 
		scale_color_discrete(name = 'Sample')+
		scale_x_continuous(breaks = c(0,200,400), labels = c('-2 Kbp','TSS','2 Kbp'))
	return(p)
}


file_list = get_profile_file_list(in_dir,pattern)

data = foreach(i = seq_along(file_list),.combine='rbind')%do%{
	fn = file_list[i]
	sample_name = names(file_list)[i]
	dt = read_profile(fn,sample_name)
	return(dt)
}

p = plot_profile(data)

fig_path = sprintf('%s/dist_to_TSS.pdf',fig_dir)
save_plot(fig_path,p,base_width=6,base_height=4)

matrix = read_matrix(matrix_fn)
row_data = matrix[[1]]
matrix = as.matrix(matrix[[2]])
weight = c(rep(0,100),rep(1,200),rep(0,100))
rank = matrix %*% weight
matrix_ordered = matrix[order(rank),]


fig_path = sprintf('%s/1346_rep1_heatmap.pdf',fig_dir)
pdf(fig_path,height=6,width=6)
heatmap3(matrix_ordered,col = topo.colors(256),Rowv=NA,Colv=NA,labRow=FALSE,labCol=FALSE,margins=c(1,1))
dev.off()