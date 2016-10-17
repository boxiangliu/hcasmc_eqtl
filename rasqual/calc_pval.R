library(data.table)


args=commandArgs(T)
in_file=args[1]
out_file=args[2]

# calculate p-value:
res=fread(in_file)
setnames(res,c('fid','rsid','chr','pos','ref','alt','freq','hwe_chisq','impute_qual','log10_qval','chisq','pi','delta','phi','overdispersion','sid','n_fsnp','n_rsnp','n_iter_null','n_iter_alt','tie','loglik_null','converge','r2_fsnp','r2_rsnp'))
res[,pval:=pchisq(q=chisq,df=1,lower.tail=F)]
res[,rank:=rank(res$pval)]


# figure: 
# pdf(out_fig)
# plot(res$pos,-log10(res$pval),col=ifelse(res$rank<=3,'red','black'),xlab='Position',ylab='-log10(P-value)')
# text(res$pos,-log10(res$pval),labels=ifelse(res$rank<=3,res$rsid,''),adj=1.1)
# dev.off()


# write metal file:
write.table(res,file=out_file,quote=F,sep='\t',row.names=F,col.names=T)
