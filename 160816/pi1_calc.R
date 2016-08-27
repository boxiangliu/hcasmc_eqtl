library(qvalue)
args=commandArgs(T)
input_file=args[1]
tissue=args[2]
pval=fread(input_file)%>%select(4)%>%unlist()
pi1=1-qvalue(pval)$pi0
writeLines(paste(tissue,pi1,sep='\t'))

