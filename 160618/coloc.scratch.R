# library:
library(coloc)
library(gtools)

# command args:
gwas_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/gwas/cad.add.160614.website.txt'
eqtl_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160530/fastqtl_nominal/fastqtl.allpairs.pc3.peer8.txt'
sdY_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160618/rpkm_sd.pc3.peer8.txt'


# read gwas data:
gwas=fread(gwas_file,header=T)%>%
	mutate(snp=paste(chr,bp_hg19,sep="_"),maf=pmin(effect_allele_freq,1-effect_allele_freq))%>%
	dplyr::select(maf,snp,beta,se_dgc,p_dgc)%>%
	rename(varbeta=se_dgc,pval=p_dgc)


# read eQTL data:
eqtl=fread(eqtl_file)%>%
	dplyr::select(gene_id=V1,snp=V2,beta=V5,varbeta=V6,pval=V4)
eqtl$snp=eqtl$snp%>%str_replace(.,'chr','')%>%str_replace(.,'_[ATGC]_[ATGC]$','')


# read expression standard deviation: 
sdY=fread(sdY_file,header=T)


# run coloc:
gene_ids=unique(eqtl$gene_id)
res=list()
for (id in gene_ids){
	eqtl2=eqtl[gene_id==id,]
	res[[id]]=coloc.abf(dataset1=list(snp=gwas$snp,beta=gwas$beta,varbeta=gwas$varbeta,type="cc",MAF=gwas$maf),
		dataset2=list(snp=eqtl2$snp,beta=eqtl2$beta,varbeta=eqtl2$varbeta,type="quant",sdY=sdY[gene_id==id,sd]))
}

# debug: 
which(names(res)%in%gene_ids)


id='ENSG00000118526.6'
eqtl2=eqtl[gene_id==id,]
res=coloc.abf(dataset1=list(snp=gwas$snp,beta=gwas$beta,varbeta=gwas$varbeta,type="cc",MAF=gwas$maf),
	dataset2=list(snp=eqtl2$snp,beta=eqtl2$beta,varbeta=eqtl2$varbeta,type="quant",sdY=sdY[gene_id==id,sd]))

which(str_detect(gwas$snp,'134214525'))
# 4020425
gwas[4020420:4020430,]


idx=which(gwas$snp%in%eqtl2$snp)
gwas2=gwas[idx,]
merged=merge(eqtl2,gwas2,by='snp')
idx=match(mixedsort(merged$snp),merged$snp)
merged=merged[idx,]
pdf()
plot(-log10(merged$pval.y),ylab='-log10(p-value)',main='TCF21')
points(-log10(merged$pval.x),col='red')
legend('topright',col=c('black','red'),legend=c('GWAS','eQTL'),pch=1)
dev.off()

which(str_detect(merged$snp,'134214525'))
# 1843
merged2=merged[1653:1893]
plot(-log10(merged2$pval.y),ylab='-log10(p-value)',main='TCF21')
points(-log10(merged2$pval.x),col='red')
legend('topright',col=c('black','red'),legend=c('GWAS','eQTL'),pch=1)
