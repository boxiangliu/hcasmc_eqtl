library(data.table)
library(cowplot)

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

color_map=get_color_map()
smr <- fread("../processed_data/finemap/smr/SMAD3_across_tissues.txt",data.table=F)
tissue <- smr$tissue
tissue <- gsub("_lite","",tissue)
tissue <- unname(tissue_id2name(tissue))
mean <- smr$b_SMR
lower <- smr$b_SMR - 1.96*smr$se_SMR
upper <- smr$b_SMR + 1.96*smr$se_SMR

df <- data.frame(tissue, mean, lower, upper)
df <- df[order(df$mean, decreasing=FALSE),]

# reverses the factor level ordering for labels after coord_flip()
df$tissue <- factor(df$tissue, levels=df$tissue)

g <- ggplot(data=df, aes(x=tissue, y=mean, ymin=lower, ymax=upper, color=tissue)) +
	geom_pointrange() + 
	geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
	scale_color_manual(values=color_map,guide='none')+
	coord_flip() +  # flip coordinates (puts labels on y axis)
	xlab("Tissue") + ylab("beta SMR (95% CI)") +
	ggtitle("SMAD3") +
	theme(plot.title = element_text(hjust = 0.8)) +
	theme_bw()  # use a white background
print(g)
