library(data.table)

complement=function(allele){
	if (allele%in%c('A','T')){
		complement=setdiff(c('A','T'),allele)
	} else {
		complement=setdiff(c('G','C'),allele)
	}
	return(complement)
}

gwas_fn='../data/gwas/howson_2017/Howson-JMM_CHD_Mixed_2017.txt'
vcf_fn='../processed_data/finemap/finemap/tmp/Howson.norm.vcf'
out_fn='../data/gwas/howson_2017/Howson-JMM_CHD_Mixed_2017.norm.txt'

gwas=fread(gwas_fn)
vcf=fread(vcf_fn)

# Append ref and alt field: 
vcf[,chrpos_b37:=sprintf('chr%s:%s',`#CHROM`,POS)]
stopifnot(vcf$chrpos_b37==gwas$chrpos_b37)
gwas[,c('ref','alt'):=list(vcf$REF,vcf$ALT)]


# Remove duplicate lines: 
gwas=gwas[!duplicated(gwas[,list(chrpos_b37,effect_allele,other_allele)]),]


# Fixe strand flips: 
temp=gwas[,setequal(c(ref,alt),c(effect_allele,other_allele)),by='chrpos_b37']
strand_flip=temp[V1==FALSE,chrpos_b37]
for (snp in strand_flip){
	gwas[chrpos_b37==snp,alt:=complement(setdiff(c(effect_allele,other_allele),complement(ref)))]
}

gwas[,snpID:=sprintf('%s:%s:%s',chrpos_b37,ref,alt)]
fwrite(gwas,out_fn,sep='\t')