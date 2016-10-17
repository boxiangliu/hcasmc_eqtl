# 
readAtac=function(command){
	atac=fread(command)
	setnames(atac,c('chrom','chromStart','chromEnd','name','score','strand','signalValue','pValue','qValue','peak'))
	message(nrow(atac),' peaks.')
	return(atac)
}


#
overlap=function(atac,gwas){
	setkey(atac,chrom,chromStart,chromEnd)
	overlaps=foverlaps(gwas,atac,nomatch=0)
	overlaps=overlaps%>%mutate(pos=i.chromStart)%>%select(-i.chromStart,-i.chromEnd)
	message(nrow(overlaps),' gwas SNPs in chromatin accessible regions.')
	return(overlaps)
}


# 
peakSummary=function(overlaps){
	overlaps_peakSummary=overlaps%>%group_by(chrom,chromStart,chromEnd)%>%summarize(n_snps=length(id))
	return(overlaps_peakSummary)
}
