library(data.table)
library(cowplot)

fig_dir = '../figures/finemap/locuszoom/gtex_comparison_smr/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

get_color_map=function(){
	color=fread('shared/tissue_color.txt')
	color[,tissue_color_hex:=max(tissue_color_hex),by=tissue]
	color_map=color$tissue_color_hex
	names(color_map)=color$tissue_site_detail
	return(color_map)
}

tissue_id2name=function(tissue_id){
	color=fread('shared/tissue_color.txt')
	color[tissue_site_detail_id=='hcasmc.rpkm',tissue_site_detail_id:='HCASMC']
	map=color$tissue_site_detail
	names(map)=color$tissue_site_detail_id
	return(map[tissue_id])
}

abbreviate_tissue_name=function(tissue_name){
	color=fread('shared/tissue_color.txt')
	color[tissue_site_detail_id=='hcasmc.rpkm',tissue_site_detail_id:='hcasmc']
	map=color$abbreviation
	names(map)=color$tissue_site_detail
	return(map[tissue_name])
}

color_map=get_color_map()

for (gene in c('TCF21','FES','SMAD3','SIPA1','PDGFRA')){
	fn = sprintf("../processed_data/finemap/smr/%s_across_tissues.txt",gene)
	smr <- fread(fn,data.table=F)
	tissue <- smr$tissue
	tissue <- gsub("_lite","",tissue)
	tissue <- unname(tissue_id2name(tissue))
	mean <- smr$b_SMR
	lower <- smr$b_SMR - 1.96*smr$se_SMR
	upper <- smr$b_SMR + 1.96*smr$se_SMR

	df <- data.frame(tissue, mean, lower, upper, stringsAsFactors = FALSE)
	df <- df[order(df$mean, decreasing=FALSE),]
	df$abbreviation=abbreviate_tissue_name(df$tissue)
	df$abbreviation <- factor(df$abbreviation, levels=df$abbreviation)

	g <- ggplot(data=df, aes(x=abbreviation, y=mean, ymin=lower, ymax=upper, color=tissue)) +
		geom_pointrange() + 
		geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
		scale_color_manual(values=color_map,guide='none')+
		ggtitle(gene) +
		coord_flip() + 
		theme(axis.text.y=element_text(color=ifelse(df$abbreviation=='HCASMC','purple','black')))+
		xlab("") + ylab("beta SMR (95% CI)")

	fig_fn = sprintf('%s/gtex_comparison_smr.%s.pdf',fig_dir,gene)
	save_plot(fig_fn,g,base_width = 6,base_height = 7.5)
}
