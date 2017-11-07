# Find a HCASMC specific peak around TCF21
library(data.table)
library(dplyr)
library(dtplyr)

args=commandArgs(T)
in_fn=args[1]
in_fn='../processed_data/hcasmc_specific_open_chromatin/intersect/intersect.bed'

read_multi_intersect=function(in_fn){
	x=fread(in_fn)
	return(x)
}


read_gencode=function(gencode_fn='../data/gtex/gencode.v19.genes.v6p.hg19.bed'){
	gencode=fread(gencode_fn,col.names=c('chr','start','end','strand','gene_id','gene_name','type'))
	return(gencode)
}



subset_to_gene=function(x,gname,window,n_tissue=1,annotation=gencode){
	window=annotation[gene_name==gname,list(chr,start=start-window,end=end+window)]
	print(window)
	y=x%>%filter(chrom==window$chr,start>=window$start,end<=window$end)%>%
	filter(HCASMC==1,num<n_tissue)

	return(y)
}

multi_intersect=fread(in_fn)
gencode=read_gencode()
tcf21=subset_to_gene(multi_intersect,'TCF21',5e5,n_tissue=3)
tcf21




find_specific_peak=function(x){
	return(0)
}

main=function(){

}