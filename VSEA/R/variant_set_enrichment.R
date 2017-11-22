# Perform VSEA
# Boxiang Liu (bliu2@stanford.edu)
# 2017-11-15

library(sinib)
overlap=function(ld_set,region){
	message('INFO - performing overlap')
	ld_set[,c('start','end'):=list(pos,pos)]

	stopifnot(all(c('chr','start','end')%in%names(region)))
	setkey(region,chr,start,end)

	# Overlap: 
	overlap=unique(foverlaps(ld_set,region[,list(chr,start,end)]))
	overlap[,c('i.start','i.end'):=NULL]

	overlap[,snp_overlap:=!is.na(start)]
	overlap[,loci_overlap:=any(snp_overlap),by='background_variant']
	overlap[,c('start','end'):=NULL]
	overlap=unique(overlap)
	stopifnot(nrow(overlap)==nrow(ld_set))
	stopifnot(overlap$snpID==ld_set$snpID)


	overlap=overlap[ld_proxy==FALSE,]
	overlap[,p:=mean(loci_overlap),by='foreground_variant']
	return(overlap)
}

calc_enrichment=function(overlap){
	# Calculate enrichment p-value:
	message('INFO - calculating enrichment p-value')
	p=overlap[,list(p=unique(p)),by='foreground_variant']
	n=as.integer(rep(1,length(p$p)))
	s=sum(overlap[snpID==foreground_variant,loci_overlap])
	psinib(s,n,p$p,lower.tail=FALSE)
}

