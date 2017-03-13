library(data.table)
library(cowplot)
library(stringr)

# Variables: 
in_file='../processed_data/hcasmc_specific_open_chromatin/intersect/intersect.quant_norm.count'
out_dir='../processed_data/hcasmc_specific_open_chromatin/specificity_index/'


# Functions:
remove_zero_rows=function(x){
	y=x[rowSums(x)!=0,]
	return(y)
}


calculate_esi=function(x,tissue){
	entropy=function(x){
		w=apply(x,1,function(y) {z=y[y!=0];-sum(z*log2(z))})
		return(w)
	}

	# toi=tissue of interest
	# ot=other tissues
	stopifnot(class(x)=='data.frame')

	toi=x[,tissue]
	ot=x[,colnames(x)!=tissue]

	print('normalizing...')
	rs=rowSums(ot)
	toin=toi/rs
	otn=ot/rs

	print('calculating specificity...')
	e=entropy(otn)
	z=e-log2(toin)

	print('removing infinity rows...')
	z=z[!is.infinite(z)]
	z=z[!is.na(z)]
	esi=1-z/max(z)

	return(esi)
}


calculate_esi_tissue_list=function(x,list_){
	esi_ls=list()
	stopifnot(all(list_%in%colnames(x)))
	for (i in list_){
		print(i)
		esi_ls[[i]]=calculate_esi(x,i)
	}
	return(esi_ls)
}


make_qqplot=function(esi_ls,plot=T){
	tmpls=list()
	for (i in names(esi_ls)){
		z=qqplot(esi_ls[[i]],esi_ls[['HCASMC']],plot.it=F)
		tmpls[[i]]=data.frame(x=z$x,y=z$y,sample=i)
	}
	qqdata=Reduce(rbind,tmpls)
	# p=ggplot(qqdata,aes(x=y,y=x,color=sample))+geom_line(alpha=ifelse(qqdata$sample=='HCASMC',1,0.2))+scale_color_manual(guide='none',breaks=tissue_color$tissue,values=tissue_color$color)+xlab('HCASMC quantiles')+ylab('GTEx quantiles')
	if (plot){
		p=ggplot(qqdata,aes(x=y,y=x,color=sample))+geom_line(alpha=ifelse(qqdata$sample=='HCASMC',1,0.2))+xlab('HCASMC quantiles')+ylab('GTEx quantiles')
		return(p)
	} else {
		return(qqdata)
	}
}


parse_id=function(id){
	parsed_id=data.frame(str_split_fixed(id,'_',4)[,1:3])
	colnames(parsed_id)=c('chrom','chromStart','chromEnd')
	return(parsed_id)
}


make_output=function(si,tissue,parsed_id){

	# Turn specificity index into a data.frame:
	si_df=data.frame(id=names(si[[tissue]]),psi=si[[i]])
	

	# Merge specificity with parsed id: 
	out=merge(parsed_id,si_df,by='id',all.x=F,all.y=T)
	

	# Rearrange columns:
	out=out[,c('chrom','chromStart','chromEnd','id','psi')]
	return(out)
}


# Create dir if necessary: 
if (dir.exists(out_dir)){dir.create(out_dir)}


# Read signal values: 
sig=fread(in_file)


# Split signal into coldata, rowdata, and matdata: 
rowdata=sig$id
matdata=as.data.frame(sig[,2:ncol(sig)])
rownames(matdata)=rowdata
coldata=colnames(matdata)


# Remove rows with all zeros: 
matdata_filt=remove_zero_rows(matdata)


# Calculate specificity index for each tissue: 
si=calculate_esi_tissue_list(matdata_filt,coldata)


# Plot the distribution of specificity index: 
p=make_qqplot(si)
qqdata=make_qqplot(si,F)
p=ggplot(qqdata,aes(x=y,y=x,color=sample))+geom_line(aes(alpha=ifelse(qqdata$sample=='HCASMC',1,0.3)))+xlab('HCASMC quantiles')+ylab('GTEx quantiles')+scale_alpha_continuous(guide='none')+scale_color_discrete(guide='none')


# Save plot:
save_plot('../figures/hcasmc_specific_open_chromatin/si_distribution_qqplot.pdf',p)


# Parse ID into chrom, chromStart, chromEnd:
parsed_id=parse_id(rowdata)
parsed_id=cbind(data.frame(id=rowdata),parsed_id)


# Output the specificity index: 
# Format: 
# chrom chromStart chromEnd id specificity_index
for (i in coldata){
	print(i)
	out=make_output(si,i,parsed_id)
	tmp_fn=tempfile()
	out_fn=sprintf("%s/%s.bed",out_dir,i)
	fwrite(out,tmp_fn,sep='\t',col.names=F)
	system(sprintf('sort -k1,1 -k2,2n %s > %s',tmp_fn,out_fn))
	file.remove(tmp_fn)
}
