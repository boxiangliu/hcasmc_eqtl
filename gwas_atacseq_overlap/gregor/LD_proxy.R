library(data.table)
library(dplyr)
library(dtplyr)
library(stringr)
library(foreach)
library(doMC)
registerDoMC(cores=11)

# Variables:
vcf_dir='../processed_data/gwas_atacseq_overlap/prepare_vcf/'
tmp_dir='../processed_data/gwas_atacseq_overlap/tmp/'
if (!dir.exists(tmp_dir)) {dir.create(tmp_dir)}
window=1e6


# Functions: 
liftOver=function(x,chain='/srv/persistent/bliu2/shared/chain_files/hg38ToHg19.over.chain.gz'){
	x[,rowid:=1:nrow(x)]
	x[,CHR_POS:=as.integer(CHR_POS)]
	bedold=x[,.(CHR_ID,CHR_POS-1,CHR_POS,rowid)]
	bedold[,CHR_ID:=paste0('chr',CHR_ID)]
	dir.create('/srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_gene/')
	bedold_file='/srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_gene/gwas_catalog.hg38.bed'
	bednew_file='/srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_gene/gwas_catalog.hg19.bed'
	fwrite(bedold,bedold_file,sep='\t',quote=F,col.names=F)
	system(paste('/srv/persistent/bliu2/tools/ucsc_tools/liftOver',bedold_file,chain,bednew_file,'/srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_gene/gwas_catalog.unmapped.bed'))
	bednew=fread(bednew_file)[,c(3,4)]
	setnames(bednew,c('CHR_POS_19','rowid'))
	y=merge(bednew,x,by='rowid')
	y[,rowid:=NULL]
	unlink('/srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_gene/',recursive=T)
	return(y)
}


rmdupvar=function(x){
	y=x[,.(CHR_ID,CHR_POS_19,`P-VALUE`)]
	y=as.data.table(y%>%group_by(CHR_ID,CHR_POS_19)%>%summarise(`P-VALUE`=min(`P-VALUE`)))
	z=merge(x,y,by=c('CHR_ID','CHR_POS_19','P-VALUE'))
	return(z)
}

select_background_variants=function(x,snpsnap){
	# INPUT: 
	# x: a list of snp IDs (chr:pos)

	background_ls=list()
	for (i in 1:length(x)){

		snpid=x[i]
		print(sprintf('INFO - %s',snpid))
		tmp=snpsnap[snpID==snpid,list(freq_bin,dist_nearest_gene_snpsnap_protein_coding,friends_ld07,gene_count)]
		freq=tmp$freq_bin
		dist=tmp$dist_nearest_gene_snpsnap_protein_coding
		ld07=tmp$friends_ld07
		count=tmp$gene_count

		step=0
		tmp=data.frame()
		while (nrow(tmp)<501){
			step=step+1
			ub=1+0.1*step
			lb=1-0.1*step
			tmp=snpsnap[(freq_bin<=freq+step)&(freq_bin>=freq-step)&(gene_count<=ub*count)&(gene_count>=lb*count)&(dist_nearest_gene_snpsnap_protein_coding<=ub*dist)&(dist_nearest_gene_snpsnap_protein_coding>=lb*dist)&(friends_ld07<=ub*ld07)&(friends_ld07>=lb*ld07),]
			print(sprintf('INFO - %s background variants selected',nrow(tmp)))
		}
		print(sprintf('INFO - tolerance: %s',step))
		set.seed(42)
		tmp=tmp[!snpID%in%x,]
		if (nrow(tmp)==500){
			background=tmp
		} else {
			background=tmp[sample(1:nrow(tmp),500),]
		}
		background_ls[[snpid]]=background
	}
	return(background_ls)
}



# Read NHGRI GWAS catalog:
gwas=fread('../data/gwas/nhgri_catalog/gwas_catalog_v1.0.1-associations_e87_r2017-01-23.tsv')


# Select variants from Nikpay et al:
cad=gwas[`FIRST AUTHOR`=='Nikpay M']
cad=liftOver(cad)
cad=rmdupvar(cad)
cad=cad[,.(CHR_ID,CHR_POS_19,`P-VALUE`,PUBMEDID,`FIRST AUTHOR`,`DISEASE/TRAIT`)]
cad[,snpID:=paste(CHR_ID,CHR_POS_19,sep=':')]


# Read SNPsnap database: 
snpsnap=fread('/srv/persistent/bliu2/shared/SNPsnap/kb1000_collection.tab',select=c(1,2,3,4,5,7,25:29))
snpsnap[,chr:=str_split_fixed(snpID,':',2)[,1]]
snpsnap=snpsnap[chr%in%c(1:22),]


# Select background variants:
background_ls=select_background_variants(cad$snpID,snpsnap)
saveRDS(background_ls,sprintf('%s/background_list.rds',tmp_dir))


# Create set of index SNPs:
container=list()
for (gwas_index in names(background_ls)){
	tmp=rbind(snpsnap[snpID==gwas_index,],background_ls[[gwas_index]])
	tmp$gwas_index=gwas_index
	container[[length(container)+1]]=tmp
}
index_set=Reduce(rbind,container)
index_set[,chr:=str_split_fixed(snpID,':',2)[,1]]
index_set[,pos:=as.integer(str_split_fixed(snpID,':',2)[,2])]
index_set[,c('window_start','window_end'):=list(pos-window,pos+window)]
fwrite(index_set,sprintf('%s/index_set.tsv',tmp_dir),sep='\t')


# Calculate LD r2:
foreach(c=1:22)%dopar%{
	print('INFO - calculating LD...')
	ld_out_prefix=sprintf('%s/ld0.7.chr%s',tmp_dir,c)
	command=sprintf('plink --vcf %s/chr%s.vcf.gz --keep-allele-order --r2 --ld-snps %s --ld-window-kb 1000 --ld-window-r2 0.7 --out %s',vcf_dir,c,paste(unique(index_set[chr==c,snpID]),collapse=','),ld_out_prefix)
	system(command)
}
command=sprintf('cat %s/ld0.7.chr*.ld | grep -v "CHR_A" > %s/ld0.7.all_chr.ld',tmp_dir,tmp_dir)
system(command)


# Get Create LD snp set:
print('INFO - creating LD snp set...')
ld=fread(sprintf('%s/ld0.7.all_chr.ld',tmp_dir),col.names=c('CHR_A','BP_A','SNP_A','CHR_B','BP_B','SNP_B','R2'))
container=list(nrow(index_set))
n=0
for (i in 1:nrow(index_set)){
	n=n+1
	index_snp=index_set[i,snpID]
	gwas_index=index_set[i,gwas_index]


	if (n%%1000==1){
		print(sprintf('INFO - %.2f%% completed',100*n/nrow(index_set)))
	}

	ld_snp=data.table(snpID=ld[SNP_A==index_snp,SNP_B],r2=ld[SNP_A==index_snp,R2])
	ld_snp$loci_index=index_snp
	ld_snp$gwas_index=gwas_index
	ld_snp[,ld_proxy:=snpID!=index_snp]
	container[[i]]=ld_snp
}
ld_set=Reduce(rbind,container)
ld_set[,chr:=paste0('chr',str_split_fixed(snpID,':',2)[,1])]
ld_set[,pos:=as.integer(str_split_fixed(snpID,':',2)[,2])]


stopifnot(length(unique(ld_set$gwas_index))==65)
stopifnot(ld_set[ld_proxy==FALSE,all(snpID==loci_index)])
stopifnot(ld_set[ld_proxy==TRUE,all(snpID!=loci_index)])
fwrite(ld_set,sprintf('%s/ld_set.tsv',tmp_dir),sep='\t')