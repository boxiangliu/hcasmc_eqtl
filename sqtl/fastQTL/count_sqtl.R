# Count the number of sQTLs:

library(data.table)
library(foreach)

args=commandArgs(T)
in_fn=args[1]
out_dir=args[2]
out_prefix=args[3]
if (!dir.exists(out_dir)){dir.create(out_dir,recursive=TRUE)}

# Functions:
fread2=function(fn,...){

	if (grepl('.gz',fn)){
		x=fread(sprintf('zcat %s',fn),...)
	} else {
		x=fread(fn,...)
	}
	return(x)
}


count_field=function(fn){
	ncol(fread2(fn,nrow=1))
}


select_top=function(x,group_by='cluster',select_by='pval'){
	setnames(x,c(group_by,select_by),c('cluster','pval'))
	
	x[,pmin:=min(pval),by='cluster']
	y=x[pval==pmin]
	stopifnot(all(y$pval==y$pmin))

	x[,pmin:=NULL]
	y[,pmin:=NULL]

	setnames(x,c('cluster','pval'),c(group_by,select_by))
	setnames(y,c('cluster','pval'),c(group_by,select_by))

	return(y)
}

main=function(in_fn,out_dir,out_prefix,fdr=c(0.001,0.01,0.05)){

	n_field=count_field(in_fn)

	message('INFO - reading input')
	if (n_field==6){

		message('INFO - input is in FastQTL nominal format')
		sqtl=fread2(in_fn,col.names=c('cluster','snp','dist','pval','beta','se'))
		sqtl=select_top(sqtl)

	} else if(n_field==16){

		message('INFO - input is in FastQTL permutation format')
		sqtl=fread2(in_fn,select=c(1,6,7,16),col.names=c('cluster','snp','dist','pval'))

	} else {

		stop('INFO - unknown format!')

	}

	sqtl[,padj:=p.adjust(pval,'fdr')]

	sig=foreach(i=fdr,.combine='rbind')%do%{
		sig=sum(unlist(sqtl[,padj<i]),na.rm=TRUE)
		data.table(fdr=i,sig=sig)
	}

	fwrite(sig,sprintf('%s/%s.sig.txt',out_dir,out_prefix),sep='\t')

}

# Main: 
main(in_fn,out_dir,out_prefix)

