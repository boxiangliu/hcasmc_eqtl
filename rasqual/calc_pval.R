library(data.table)


args=commandArgs(T)
in_file=args[1]
out_file=args[2]


# Skip existing files: 
if (file.exists(out_file)){
	stop(paste(out_file, " exists. Skipping..."))
}

# calculate p-value:
res=fread(in_file)
setnames(res,c('fid','rsid','chr','pos','ref','alt','freq','hwe_chisq','impute_qual','log10_qval','chisq','pi','delta','phi','overdispersion','sid','n_fsnp','n_rsnp','n_iter_null','n_iter_alt','tie','loglik_null','converge','r2_fsnp','r2_rsnp'))
res[,pval:=pchisq(q=chisq,df=1,lower.tail=F)]
res[,rank:=rank(res$pval)]


# output file:
write.table(res,file=out_file,quote=F,sep='\t',row.names=F,col.names=T)
