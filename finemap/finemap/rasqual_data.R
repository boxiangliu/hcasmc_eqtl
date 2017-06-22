library(data.table)
library(stringr)
library(cowplot)
library(foreach)
library(doMC)
registerDoMC(cores=20)


gwas_fn='../data/gwas/howson_2017/Howson-JMM_CHD_Mixed_2017.norm.in1kgp3.txt'
gencode_fn='/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/gencode.v19.genes.v6p.hg19.bed'
eqtl_dir='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/rasqual/output_pval/'
ld_block_fn='/users/bliu2/tools/LDetect/ldetect-data/EUR/fourier_ls-all.bed'
vcf_dir='/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/asvcf_sid/'
out_dir='../processed_data/finemap/finemap/rasqual_data/input/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}


# Functions: 
calc_zscore=function(pval,pi){
	((pi>=0.5)*2-1)*abs(qnorm(pval/2))
}


# Main: 
gwas=fread(gwas_fn,header=TRUE)

temp=gwas[,str_split_fixed(chrpos_b37,':',2)]
gwas[,c('chr','pos'):=list(temp[,1],temp[,2])]
gwas[,c('start','end'):=as.integer(pos)]


ld_block=fread(ld_block_fn,header=TRUE,
	col.names=c('chr','start','end'))
ld_block[,start:=start+1]
ld_block[,lociID:=sprintf('%s_%s_%s',chr,start,end)]

setkey(gwas,chr,start,end)
setkey(ld_block,chr,start,end)
overlap=foverlaps(gwas[,list(chr,start,end,rsid,snpID,
	ref,alt,p,n,chrpos_b37)],
	ld_block,nomatch=0)

overlap[,pos:=i.start]
overlap[,c('i.start','i.end'):=NULL]

min_p=overlap[,list(min_p=min(p)),by='lociID']
sig_loci=min_p[min_p<1e-6,lociID]
sig_loci=data.table(sig_loci,str_split_fixed(sig_loci,'_',3))
setnames(sig_loci,c('lociID','chr','start','end'))
sig_loci[,c('start','end'):=list(as.integer(start),as.integer(end))]

gencode=fread(gencode_fn)
setnames(gencode,c('chr','start','end',
	'strand','gene_id','gene_name','type'))

window_size=1e6
gencode[,window_start:=ifelse(strand=='+',start-window_size,end-window_size)]
gencode[,window_end:=ifelse(strand=='+',start+window_size,end+window_size)]

setkey(gencode,chr,window_start,window_end)
setkey(sig_loci,chr,start,end)
overlap=foverlaps(sig_loci,gencode,nomatch=0)
overlap=overlap[type%in%c('protein_coding','lincRNA')]

eqtl_threshold=1e-6
container=list()
n=0
for (i in 1:nrow(overlap)){
	message(i)
	chr=overlap[i,chr]
	gene_id=overlap[i,gene_id]
	gene_name=overlap[i,gene_name]

	eqtl_fn=sprintf('%s/%s/%s_%s.pval.txt',eqtl_dir,chr,gene_id,gene_name)

	eqtl=fread(eqtl_fn)
	if (min(eqtl$pval)<eqtl_threshold){
		n=n+1
		container[[n]]=eqtl
	}
}

eqtl=Reduce(rbind,container)
eqtl=eqtl[!(duplicated(eqtl[,list(chr,pos,fid,rsid)]) | 
	duplicated(eqtl[,list(chr,pos,fid,rsid)], fromLast=TRUE)),]


chrpos_fn=sprintf('%s/eqtl.chrpos.txt',out_dir)
fwrite(unique(eqtl[,list(chr,pos)]),chrpos_fn,col.names=F,sep='\t')
system('bash finemap/finemap/eqtl_snps.sh')

vcf=fread('cut -f1-5,8 ../processed_data/finemap/finemap/rasqual_data/eqtl.vcf',skip=45)
setnames(vcf,c('chr','pos','rsid','ref','alt','info'))
vcf=unique(vcf)


vcf=vcf[!str_detect(alt,','),]

vcf[,freq:=as.numeric(str_extract(info,'(?<=AF=)(.+?)(?=;)'))]


merged=merge(eqtl,vcf,by=c('chr','pos','rsid'))
merged[,freq_diff:=abs(freq.x-freq.y)]
merged2=merged[(ref.x==ref.y&alt.x==alt.y)|freq_diff<=0.05,list(chr,pos,rsid,fid,ref=ref.y,alt=alt.y,pval,pi)]
merged2[,zscore:=calc_zscore(pval,pi)]


merged2[,c('start','end'):=list(as.integer(pos),as.integer(pos))]
setkey(merged2,chr,start,end)
overlap=foverlaps(merged2,sig_loci,nomatch=0)
overlap[,pos:=i.start]
overlap[,c('i.start','i.end'):=NULL]
overlap[,snpID:=paste(chr,pos,ref,alt,'b37',sep='_')]
overlap[,snpID:=str_replace(snpID,'chr','')]

iterator=unique(overlap[,list(lociID,fid)])

for (i in 1:nrow(iterator)){
	print(i)
	id=iterator[i,lociID]
	gene_id=iterator[i,fid]
	output=overlap[lociID==id&fid==gene_id,list(snpID,zscore)]
	out_fn=sprintf('%s/%s.%s.z',out_dir,id,gene_id)
	write.table(output,out_fn,col.names=F,row.names=F,quote=F,sep=' ')
}

print('INFO - calculating LD...')
vcf_fn='../processed_data/finemap/finemap/rasqual_data/eqtl.id.vcf'
foreach(i = 1:nrow(iterator))%dopar%{

	print(sprintf('INFO - %s',id))
	id=iterator[i,lociID]
	gene_id=iterator[i,fid]


	snps=overlap[lociID==id&fid==gene_id,list(lociID,snpID)]
	snps[,snpID:=str_replace(snpID,'chr','')]


	snp_fn=sprintf('%s/%s.%s.snpid',out_dir,id,gene_id)
	fwrite(snps[,2], snp_fn, col.names=F)

	out_prefix=sprintf('%s/%s.%s',out_dir,id,gene_id)
	command=sprintf('plink --vcf %s --extract %s --keep-allele-order --r square --out %s',vcf_fn,snp_fn,out_prefix)
	system(command)
	
	file.remove(paste0(out_prefix,'.nosex'),paste0(out_prefix,'.log'))
}


# Sanity checks: 
stopifnot(length(list.files(out_dir,'ld'))==length(list.files(out_dir,'z')))
for (i in 1:nrow(iterator)){
	id=iterator[i,lociID]
	gene_id=iterator[i,fid]
	n1=as.integer(str_match(system(sprintf('wc -l %s/%s.%s.ld',out_dir,id,gene_id),intern=TRUE),'[0-9]+'))
	n2=as.integer(str_match(system(sprintf('wc -l %s/%s.%s.z',out_dir,id,gene_id),intern=TRUE),'[0-9]+'))
	stopifnot(n1==n2)
}
