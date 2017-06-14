library(data.table)
library(stringr)

gwas_fn='../data/gwas/howson_2017/Howson-JMM_CHD_Mixed_2017.txt'
ld_block_fn='/users/bliu2/tools/LDetect/ldetect-data/EUR/fourier_ls-all.bed'
out_dir='../processed_data/finemap/enloc/'

if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

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
overlap=foverlaps(gwas[,list(chr,start,end,chrpos_b37,zscore)],
	ld_block,nomatch=0)

output=overlap[,list(chrpos_b37,lociID,zscore)]
out_fn=gzfile(sprintf('%s/howson.zscore.gz',out_dir),'w')
write.table(output,out_fn,col.names=F,row.names=F,quote=F,sep='\t')





