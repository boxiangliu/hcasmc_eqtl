library(data.table)
library(stringr)
library(cowplot)
library(foreach)
library(doMC)
registerDoMC(cores=20)


gwas_fn='../data/gwas/howson_2017/Howson-JMM_CHD_Mixed_2017.txt'
ld_block_fn='/users/bliu2/tools/LDetect/ldetect-data/EUR/fourier_ls-all.bed'
vcf_dir='../processed_data/gwas_atacseq_overlap/prepare_vcf/'
out_dir='../processed_data/finemap/finemap/'
fig_dir='../figures/finemap/finemap/'
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
overlap=foverlaps(gwas[,list(chr,start,end,rsid,
	effect_allele,other_allele,n,chrpos_b37,zscore)],
	ld_block,nomatch=0)
overlap[,pos:=i.start]
overlap[,c('i.start','i.end'):=NULL]

N=overlap[,list(.N),by='lociID']
p=ggplot(N,aes(N))+geom_histogram()
fig_fn=sprintf('%s/howson.snp_per_loci.pdf',fig_dir)
save_plot(fig_fn,p)

for (id in unique(overlap$lociID)){
	print(sprintf('INFO - %s',id))
	output=overlap[lociID==id,list(chrpos_b37,zscore)]
	out_fn=sprintf('%s/%s.z',out_dir,id)
	write.table(output,out_fn,col.names=F,row.names=F,quote=F,sep='\t')
}


# Calculate LD r2:
print('INFO - calculating LD...')
foreach(id=unique(overlap$lociID))%dopar%{
	
	print(sprintf('INFO - %s',id))
	
	snps=overlap[lociID==id,list(chr=str_replace(chr,'chr',''),
		pos,rsid=str_replace(chrpos_b37,'chr',''),
		effect_allele,other_allele)]
	
	chr=unique(unlist(snps$chr))
	command=sprintf("bcftools query -r %s -f '%%CHROM  %%POS  %%REF  %%ALT\n' %s/chr%s.vcf.gz",
		paste(snps$rsid,collapse=','),vcf_dir,chr)
	res=system(command,intern=TRUE)

	res=data.table(str_split_fixed(res,'  ',4))
	setnames(res,c('chr','pos','ref','alt'))
	res[,pos:=as.integer(pos)]

	stopifnot(nrow(res)!=nrow(snps))

	merged=merge(res,snps,by=c('chr','pos'))
	merged[,same:=setequal(c(ref,alt),c(effect_allele,other_allele)),
		by=c('chr','pos','ref','alt')]
	
	stopifnot(all(merged$same))

	out_prefix=sprintf('%s/%s',out_dir,id)
	command=sprintf("bcftools view -r %s %s/chr%s.vcf.gz -Oz > %s/%s.vcf.gz",paste(snps$rsid,collapse=','),vcf_dir,chr,out_dir,id)
	system(command)

	command=sprintf('plink --vcf %s/%s.vcf.gz --keep-allele-order --r square --out %s',out_dir,id,out_prefix)
	system(command)
	file.remove(paste0(out_prefix,'.nosex'),paste0(out_prefix,'.log'))
}




