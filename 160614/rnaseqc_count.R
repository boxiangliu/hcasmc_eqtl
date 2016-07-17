#!/usr/bin/env Rscript
# boxiang liu
# durga
# extract read counts from rnaseqc output 

# library:
library(R.utils)


# command args:
args=commandArgs(T,T)
sample_dirs_file=args$input


# read input list:
sample_dirs=scan(sample_dirs_file,what=character())
message(length(sample_dirs),' samples.')


# read samples: 
for (sample_dir in sample_dirs){
	message(sample_dir)

	# go into sample directory:
	cwd=getwd()
	setwd(paste0(sample_dir,'/report'))
	

	# read count and rpkm:
	sample=basename(sample_dir)
	rpkm=fread('genes.rpkm.gct',header=T)
	count=fread(sprintf('%s/%s.metrics.tmp.txt.intronReport.txt',sample,sample),header=T)
	

	# merge count and rpkm: 
	count=count[,.(Transcript,Exon_Reads)]
	merged=merge(rpkm,count,by.x='Name',by.y='Transcript',all=T)
	stopifnot(nrow(rpkm)==nrow(merged))

	# prepare to write count:
	output=merged[,.(Name,Description,Exon_Reads)]
	output$Exon_Reads[is.na(output$Exon_Reads)]=0
	idx=match(rpkm$Name,output$Name)
	output=output[idx,]
	setnames(output,'Exon_Reads',sample)


	# write output header: 
	writeLines(sprintf('#1.2\n%s\t%s',nrow(output),1),'genes.reads.gct')
	write.table(output,'genes.reads.gct',append=T,quote=F,sep='\t',row.names=F,col.names=T)


	# go to original directory:
	setwd(cwd)
} 
