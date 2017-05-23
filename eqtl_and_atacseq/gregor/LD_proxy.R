library(data.table)
library(dplyr)
library(dtplyr)
library(stringr)
library(foreach)
library(doMC)
registerDoMC(cores=11)

# Variables:
in_fn='../processed_data/compare_rasqual_and_fastqtl/top_snp_per_gene/rasqual.txt'
tmp_dir='../processed_data/eqtl_and_atacseq/tmp/'
if (!dir.exists(tmp_dir)) {dir.create(tmp_dir)}


# Functions:
select_background_variants=function(x,snpsnap){
	# INPUT: 
	# x: a list of snp IDs (chr:pos)


	print(sprintf('INFO - sorting snpsnap by snpID.'))
	setkey(snpsnap,snpID)
	snps=snpsnap[CJ(x),nomatch=0L]

	print(sprintf('INFO - sorting snpsnap by attributes.'))
	setkey(snpsnap,freq_bin,gene_count,friends_ld07)

	background_ls=list()
	for (i in 1:nrow(snps)) {
		snpid=snps[i,snpID]
		freq=snps[i,freq_bin]
		dist=snps[i,dist_nearest_gene_snpsnap_protein_coding]
		ld07=snps[i,friends_ld07]
		count=snps[i,gene_count]
		print(sprintf('INFO - %s',snpid))


		step=0
		tmp=data.frame()
		while ( (nrow(tmp)<501) & (step<5)){
			step=step+1
			flb=freq-step
			fub=freq+step

			ub=1+0.1*step
			lb=1-0.1*step

			clb=round(count*lb)
			cub=round(count*ub)


			llb=round(ld07*lb)
			lub=round(ld07*ub)


			tmp=snpsnap[CJ(flb:fub,clb:cub,llb:lub),nomatch=0L]


			dlb=round(dist*lb)
			dub=round(dist*ub)


			tmp=tmp[dist_nearest_gene_snpsnap_protein_coding%in%dlb:dub,nomatch=0L]
			tmp=tmp[!snpID%in%x,]
		}
		print(sprintf('INFO - tolerance: %s',step))
		print(sprintf('INFO - %s background variants selected',nrow(tmp)))

		set.seed(42)
		if (nrow(tmp)>500){tmp=tmp[sample(1:nrow(tmp),500),]}
		if (nrow(tmp)>0){background_ls[[i]]=data.table(tmp,query=snpid)}
		rm(tmp)

	}
	return(background_ls)
}


# Read NHGRI GWAS catalog:
ras=fread(in_fn)
ras[,chr:=str_replace(chr,'chr','')]

# Read SNPsnap database: 
snpsnap=fread('/srv/persistent/bliu2/shared/SNPsnap/kb1000_collection.tab',select=c(1,2,3,4,5,7,27))
tmp=as.data.table(snpsnap[,str_split_fixed(snpID,':',2)])
setnames(tmp,c('chr','pos'))
snpsnap=cbind(snpsnap,tmp)
snpsnap=snpsnap[chr%in%c(1:22),]
snpsnap[,pos:=as.integer(pos)]
snpsnap[,rsID:=ifelse(str_detect(rsID,'rs'),rsID,'.')]
rm(tmp)


# Group eQTLs based on MAF, LD buddies, etc: 
merged=merge(ras[,list(rsid,chr,pos,pval)],snpsnap,by.x=c('rsid','chr','pos'),by.y=c('rsID','chr','pos'))


# Select background variants:
background_ls=select_background_variants(merged$snpID,snpsnap)
stopifnot(length(background_ls)==length(merged$snpID))
saveRDS(background_ls,sprintf('%s/background_list.rds',tmp_dir))
background=Reduce(rbind,background_ls)
rm(background_ls)


# Select the top 1000 eGenes:
set.seed(42)
merged[,rank:=rank(pval,ties='random')]
top=merged[rank<=1000,]
background_ls=select_background_variants(top$snpID,snpsnap)
background=Reduce(rbind,background_ls)
setnames(background,'query','eqtl')
tmp=snpsnap[snpID%in%background[,unique(eqtl)],]
tmp[,eqtl:=snpID]
background=rbind(background,tmp)
fwrite(background,sprintf('%s/background_set.top1000eqtl.tsv',tmp_dir),sep='\t')

