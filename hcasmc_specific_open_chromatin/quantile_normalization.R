library(data.table)
library(cowplot)
library(preprocessCore)

# Functions:
plot_sig_vs_rank=function(sig,downsample=1e4){
	sig_down=sig[sample(1:nrow(sig),downsample),]
	sig_long=melt(sig_down,id.vars='id',value.name='signal',variable.name='epigenome')
	sig_long[,rank:=rank(-signal),by='epigenome']

	p=ggplot(sig_long,aes(rank,log10(signal)))+geom_line(alpha=ifelse(sig_long$epigenome=='HCASMC',1,0.1))

	return(p)
}


# Constants: 
in_file='../processed_data/hcasmc_specific_open_chromatin/intersect/intersect.count'
out_file='../processed_data/hcasmc_specific_open_chromatin/intersect/intersect.quant_norm.count'
out_fig='../figures/hcasmc_specific_open_chromatin/signal_versus_rank.pdf'

# Read in peak signal: 
sig=fread(in_file,header=T)


# Plot peak signal versus rank:
p1=plot_sig_vs_rank(sig)


# Make rownames the ID column:
rowdata=sig$id
coldata=colnames(sig)[2:ncol(sig)]
sig_mat=as.matrix(sig[,2:ncol(sig)])

rownames(sig_mat)=rowdata
colnames(sig_mat)=coldata


# Perform quantile normalization: 
sig_mat_norm=normalize.quantiles(sig_mat)


# Save text:  
sig_norm=as.data.table(sig_mat_norm)
colnames(sig_norm)=coldata
sig_norm=cbind(data.table(id=rowdata),sig_norm)

fwrite(sig_norm,out_file,row.names=F,col.names=T,sep='\t')


# Plot signal versus rank after quantile norm.: 
p2=plot_sig_vs_rank(sig_norm)


# Save figures: 
pdf(out_fig)
print(p1)
print(p2)
dev.off()