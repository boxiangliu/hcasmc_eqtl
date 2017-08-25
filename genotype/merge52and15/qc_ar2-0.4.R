library(foreach)
library(doMC)
library(data.table)
library(cowplot)
registerDoMC(23)

wd='/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/'
fig_dir='../figures/genotype/merge52and15_ar2-0.4/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

res=foreach(i = c(seq(22),'X'))%dopar% {
	cmd=sprintf('bcftools view -H %s/67total_ar2-0.4/chr%s.15.vcf.gz | wc -l',wd,i)
	nline_array=as.integer(system(cmd,intern=TRUE))

	cmd=sprintf('bcftools view -H %s/67total_ar2-0.4/chr%s.52.vcf.gz | wc -l',wd,i)
	nline_wgs=as.integer(system(cmd,intern=TRUE))

	cmd=sprintf('bcftools view -H %s/67total_ar2-0.4/chr%s.vcf.comm.gz | wc -l',wd,i)
	nline_intersected=as.integer(system(cmd,intern=TRUE))

	data.frame(chr=i,array=nline_array,wgs=nline_wgs,intersected=nline_intersected)
}


line_count=Reduce(rbind,res)
setDT(line_count)
line_count[,chr:=c(seq(22),'X')]
line_count[,chr:=factor(chr,level=c(seq(22),'X'))]

line_count_melt=melt(line_count,id.var='chr',variable.name='type',value.name='n_lines')
line_count_melt[,type:=factor(type,level=c('wgs','array','intersected'))]
p1=ggplot(line_count_melt,aes(x=chr,y=n_lines,fill=type))+geom_bar(stat='identity',position='dodge')+ylab('Number of records')


line_count[,pct_wgs:=intersected/wgs*100]
p2=ggplot(line_count,aes(x=chr,y=pct_wgs,label=sprintf('%.2f%%',pct_wgs)))+geom_bar(stat='identity',position='dodge')+geom_text(angle=45,nudge_y=5)+ylab('Intersected/WGS (%)')

line_count[,pct_array:=intersected/array*100]
p3=ggplot(line_count,aes(x=chr,y=pct_array,label=sprintf('%.2f%%',pct_array)))+geom_bar(stat='identity',position='dodge')+geom_text(angle=45,nudge_y=5)+ylab('Intersected/Array (%)')


res2=foreach(i = c(seq(22),'X'))%dopar% {
	cmd=sprintf('bcftools view -H --min-af 0.05 %s/67total_ar2-0.4/chr%s.vcf.comm.gz | wc -l',wd,i)
	nline_intersected_maf=as.integer(system(cmd,intern=TRUE))

	cmd=sprintf('bcftools view -H --min-af 0.05 %s/67total_ar2-0.4/chr%s.vcf.gz | wc -l',wd,i)
	nline_pooled_maf=as.integer(system(cmd,intern=TRUE))

	data.frame(chr=i,intersected=nline_intersected_maf,pooled=nline_pooled_maf)
}

line_count2=Reduce(rbind,res2)
setDT(line_count2)
line_count2[,pct:=intersected/pooled*100]
p4=ggplot(line_count2,aes(x=chr,y=pct,label=sprintf('%.2f%%',pct)))+geom_bar(stat='identity',position='dodge')+geom_text(angle=45,nudge_y=5)


pdf(sprintf('%s/line_count.pdf',fig_dir),height=6,width=10)
p1;p2;p3;p4
dev.off()


