#!/usr/bin/env Rscript
# bosh liu
# durga
# normalize PSI

# library:
library(foreach)
library(doMC)
library(cowplot)
library(preprocessCore)
source('utils.R')
registerDoMC(cores=26)

# command args:
args=commandArgs(T,T)
junc_file=args$leafcutter
figure_file=args$figure
output_file=args$output

# junc_file='/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/leafcutter/leafcutter_perind.counts.gz'
# figure_file='/srv/persistent/bliu2/HCASMC_eQTL/figures/160627/intron_missingness_rate.pdf'
# output_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160627/leafcutter.norm.tsv'

# read count file:
if (str_detect(junc_file,'.gz')){
	junc=fread(sprintf('zcat %s',junc_file),header=T)
} else {
	junc=fread(junc_file,header=T)
}


# convert fractions into doubles:
junc2=foreach(j=2:ncol(junc),.combine=cbind) %dopar% {
	temp=unlist(junc[,j,with=F])
	names(temp)=junc[,chrom]
	sapply(temp,function(x) eval(parse(text=x)))
}
colnames(junc2)=colnames(junc)[-1]


# arrange columns alphabetically:
junc2=junc2[,sort(colnames(junc2))]


# sample missingness rate: 
miss_rate=rowMeans(t(is.na(junc2)))
miss_rate=data.frame(sample=names(miss_rate),rate=miss_rate,stringsAsFactors=F)
p=ggplot(miss_rate,aes(x=sample,y=rate))+geom_point()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+background_grid(major='x')
save_plot(figure_file,p,base_width=7,base_height=7)


# mean imputation:
junc3=junc2
for (i in 1:nrow(junc3)){
	mean_value=mean(junc3[i,],na.rm=T)
	junc3[i,is.na(junc3[i,])]=mean_value
}


# quantile normalize each sample: 
junc3=normalize.quantiles(junc3,copy=T)


# inverse rank normalize each gene:
num_row=nrow(junc3)
junc4=junc3
for (i_row in seq(num_row)){
	x=junc3[i_row,]
	x_norm=getInverseNormal(x)
	junc4[i_row,]=x_norm
}


# write output: 
colnames(junc4)=colnames(junc2)
rownames(junc4)=rownames(junc2)
output=data.table(Name=rownames(junc4), junc4)
write.table(output, file=output_file, quote=F, sep='\t', row.names=F, col.names=T)

