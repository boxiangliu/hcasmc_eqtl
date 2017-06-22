library(data.table)
library(stringr)
library(cowplot)
library(foreach)
library(doMC)
registerDoMC(cores=20)


gwas_fn='../data/gwas/howson_2017/Howson-JMM_CHD_Mixed_2017.norm.in1kgp3.txt'
ld_block_fn='/users/bliu2/tools/LDetect/ldetect-data/EUR/fourier_ls-all.bed'
vcf_fn='../processed_data/finemap/finemap/preprocess/all.final.vcf'
out_dir='../processed_data/finemap/finemap/howson_data/'
fig_dir='../figures/finemap/finemap/howson_data/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

gwas=fread(gwas_fn,header=TRUE)

temp=gwas[,str_split_fixed(chrpos_b37,':',2)]
gwas[,c('chr','pos'):=list(temp[,1],temp[,2])]
gwas[,c('start','end'):=as.integer(pos)]

gwas[,zscore:=beta/se]

ld_block=fread(ld_block_fn,header=TRUE,col.names=c('chr','start','end'))
ld_block[,start:=start+1]
ld_block[,lociID:=sprintf('%s_%s_%s',chr,start,end)]

setkey(gwas,chr,start,end)
setkey(ld_block,chr,start,end)
overlap=foverlaps(gwas[,list(chr,start,end,rsid,snpID,
	ref,alt,n,chrpos_b37,zscore)],
	ld_block,nomatch=0)

overlap[,pos:=i.start]
overlap[,c('i.start','i.end'):=NULL]


N=overlap[,list(.N),by='lociID']
p=ggplot(N,aes(N))+geom_histogram()
fig_fn=sprintf('%s/howson.snp_per_loci.pdf',fig_dir)
save_plot(fig_fn,p)

# Output number of individuals:
ns=list()
for (id in unique(overlap$lociID)){
	print(sprintf('INFO - %s',id))
	out_fn=sprintf('%s/%s.n',out_dir,id)
	n=overlap[lociID==id,round(mean(n))]
	write.table(n,out_fn,quote=F,row.names=F,col.names=F)
}


# Calculate LD r2:
print('INFO - calculating LD...')
foreach(id=unique(overlap$lociID))%dopar%{

	print(sprintf('INFO - %s',id))
	
	snps=overlap[lociID==id,list(lociID,snpID)]
	snps[,snpID:=str_replace(snpID,'chr','')]
	
	snp_fn=sprintf('%s/%s.snpid',out_dir,id)
	fwrite(snps[,2], snp_fn, col.names=F)

	out_prefix=sprintf('%s/%s',out_dir,id)
	command=sprintf('plink --vcf %s --extract %s --keep-allele-order --r square --out %s',vcf_fn,snp_fn,out_prefix)
	system(command)
	
	file.remove(paste0(out_prefix,'.nosex'),paste0(out_prefix,'.log'))
}


# Output z-score:
overlap[,snpID:=str_replace_all(snpID,':','_')]
overlap[,snpID:=str_replace_all(snpID,'chr','')]
overlap[,snpID:=paste0(snpID,'_b37')]
for (id in unique(overlap$lociID)){
	print(sprintf('INFO - %s',id))
	output=overlap[lociID==id,list(snpID,zscore)]
	out_fn=sprintf('%s/%s.z',out_dir,id)
	write.table(output,out_fn,col.names=F,row.names=F,quote=F,sep=' ')
}

# Sanity checks: 
stopifnot(length(list.files(out_dir,'ld'))==length(list.files(out_dir,'z')))
for (id in unique(overlap$lociID)){
	n1=as.integer(str_match(system(sprintf('wc -l %s/%s.ld',out_dir,id),intern=TRUE),'[0-9]+'))
	n2=as.integer(str_match(system(sprintf('wc -l %s/%s.z',out_dir,id),intern=TRUE),'[0-9]+'))
	stopifnot(n1==n2)
}





