library(data.table)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(cowplot)

genome = BSgenome.Hsapiens.UCSC.hg19
finemap_fn = '../processed_data/finemap/finemap/motif_match/Finemapping.tsv'
out_dir='../processed_data/finemap/finemap/motif_match/'
fig_dir='../figures/finemap/finemap/motif_match/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

# Functions: 
str_replace_by_pos=function(x,from,to,replacement){
	y=paste0(substr(x,1,from-1),replacement,substr(x,to+1,nchar(x)))
	return(y)
}

finemap = fread(finemap_fn)
finemap[,c('chr','pos','ref','alt'):=list(
	str_split_fixed(snpid,':',4)[,1],
	as.integer(str_split_fixed(snpid,':',4)[,2]),
	str_split_fixed(snpid,':',4)[,3],
	str_split_fixed(snpid,':',4)[,4])]


out_fn=sprintf('%s/flanking_seq.fa',out_dir)
if (file.exists(out_fn)) {file.remove(out_fn)}
window=20
for (i in 1:nrow(finemap)){
	message(i)
	id=finemap[i,id]
	gene=finemap[i,gene]
	rsid=finemap[i,rsid]
	snpid=finemap[i,snpid]

	ref=finemap[i,ref]
	alt=finemap[i,alt]

	chr=finemap[i,chr]
	pos=finemap[i,pos]
	
	ref_seq=as.character(getSeq(genome,chr,start=pos-window,
		end=pos+window))
	
	if (nchar(ref)==1 & nchar(alt)==1){
		alt_seq=str_replace_by_pos(ref_seq,21,21,alt)
	} else if (nchar(ref)==1 & nchar(alt)>1){
		insert_size=nchar(alt)-1
		alt_seq=str_replace_by_pos(ref_seq,21,21,alt)
	} else if (nchar(ref)>1 & nchar(alt)==1){
		deletion_size=nchar(ref)-1
		alt_seq=str_replace_by_pos(ref_seq,21,21+deletion_size,alt)
	} else {
		stop('variant ', chr, ':', pos, ':', ref, ':', alt, ' is not recognized!')
	}

	write(sprintf('>%s_ref;%s;%s;%s\n%s',id,gene,rsid,snpid,ref_seq),file=out_fn,append=TRUE)
	write(sprintf('>%s_alt;%s;%s;%s\n%s',id,gene,rsid,snpid,alt_seq),file=out_fn,append=TRUE)
}


cmd=sprintf('/srv/persistent/bliu2/tools/meme_4.12.0/bin/fimo --oc %s/fimo/ --thresh %.02f %s %s',
	out_dir,1.0,'../data/jaspar/nonredundant/pfm_vertebrates.meme',out_fn)

system(cmd)

fimo=fread(sprintf('%s/fimo/fimo.txt',out_dir),header=TRUE)
setnames(fimo,'# motif_id','motif_id')
temp=fimo[,str_split_fixed(sequence_name,';',4)]
fimo[,c('id','gene','rsid','snpid'):=
	list(temp[,1],temp[,2],temp[,3],temp[,4])]
fimo[,logp:=-log10(`p-value`)]

temp=fimo[,str_split_fixed(id,'_',2)]
fimo[,c('id','variant'):=list(temp[,1],temp[,2])]

fimo_cast=dcast(fimo,motif_id+motif_alt_id+start+stop+strand+id+gene+rsid+snpid~variant,value.var=c('p-value','matched_sequence'))
fimo_cast[,c('log10p','log2fc'):=list(-log10(min(`p-value_alt`,`p-value_ref`,na.rm=TRUE)),-log2(`p-value_alt`/`p-value_ref`)),by=c('motif_id','start','gene','rsid','snpid')]

p=ggplot(fimo_cast,aes(log2fc,log10p))+
	geom_point(alpha=0.5,color=ifelse(fimo_cast[,abs(log2fc)>1&log10p>4],'blue','black'))+
	xlab('-log2(Fold change)')+ylab('-log10(P-value)')+
	geom_hline(yintercept=4,linetype='dashed',color='red')+
	geom_vline(xintercept=c(-1,1),linetype='dashed',color='red')

save_plot(sprintf('%s/volcano.pdf',fig_dir),p)

temp=fimo_cast[abs(log2fc)>1&log10p>4,]
setorder(temp,gene,rsid)
fwrite(temp,sprintf('%s/top_hits.txt',out_dir),sep='\t')