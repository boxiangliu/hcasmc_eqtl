#!/usr/bin/env Rscript
# bosh liu
# durga
# 16/06/26
# plot verifyBamID "FREEMIX" column

# library:
library(cowplot)
library(data.table)


# paths: 
input_file='../processed_data/genotype/quality_control/detect_WGS_contamination/verifyBAMID.combined.tsv'
figure_path='../figures/genotype/quality_control/detect_WGS_contamination/'
old_to_new_fn = '../data/joint3/orig/old_to_new_sample_name.txt'

# Function:
update_sample_names = function(sample_names,old_to_new_fn){
	old_to_new = fread(old_to_new_fn)
	setnames(old_to_new,c('old','new'))
	old_to_new = old_to_new[old != '313605']
	sample_names[sample_names %in% old_to_new$old] = old_to_new$new
	return(sample_names)
}

# read input:
input=fread(input_file,header=F)

# setnames:
setnames(input, c('sample','contamination'))

# Update sample names:
input$sample = update_sample_names(input$sample,old_to_new_fn)


# Order by contamination porportion:
setorder(input,contamination)
input[,sample := factor(sample,sample)]


# plot:
p=ggplot(input,aes(x=sample,y=contamination))+
	geom_point()+
	geom_text(aes(label=ifelse(contamination>0.1,contamination,"")),angle=0,nudge_y=0.03)+
	background_grid(major = "x")+
	ylab('Contamination Proportion')+
	xlab('')+
	geom_hline(yintercept = 0.1, color = 'red', linetype = 'dashed') +
	coord_flip()
save_plot(paste0(figure_path,'contamination_proportion.pdf'),p,base_width=8,base_height=8)