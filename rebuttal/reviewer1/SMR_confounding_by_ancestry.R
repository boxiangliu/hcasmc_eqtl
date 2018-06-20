library(data.table)
source('rebuttal/utils.R')
library(stringr)
library(dplyr)

vcf_dir = '../data/joint3/asvcf/'
smr_fn = '../processed_data/finemap/smr/UKBB_HCASMC_eqtl_smr_results5e-05.out.smr'
genotype_pc_fn = '../processed_data/genotype/genotype_pc/genotype_pcs.52samples.tsv'
rasqual_dir = '../processed_data/rasqual/output_pval/' 

read_vcf = function(vcf_fn){
	return(vcf)
}


read_rasqual = function(rasqual_fn){
	rasqual = fread(rasqual_fn)
	setorder(rasqual,pval)
	x = rasqual[1,list(fid,rsid,chr,pos,ref,alt)]
	return(x)
}

read_smr=function(smr_fn){
	fread(smr_fn)[,list(chrom=ProbeChr,gene_name=Gene,pos=Probe_bp,y=log10(p_SMR),gwas_logp=-log10(p_GWAS),eqtl_logp=-log10(p_eQTL),method='SMR')]
}

read_genotype_pc = function(genotype_pc_fn){
	genotype_pc = read.table(genotype_pc_fn,check.names=FALSE)
	x = data.frame(t(genotype_pc[1:3,]))
	x$sample = rownames(x)
	setDT(x)
	return(x)
}

smr = read_smr(smr_fn)
covariate = read_known_covariate(covariate_fn)

gene_name = smr[y< log10(5e-5),gene_name]

fn = list.files(rasqual_dir,pattern = gene_name[1],recursive=TRUE,full.names=TRUE)
rasqual = read_rasqual(fn)
vcf_fn = sprintf('%s/phased_and_imputed.%s.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz',vcf_dir,rasqual$chr)
command = sprintf('bcftools query -r %s:%s-%s -H -f "%%ID[\t%%DS]\n" %s',rasqual$chr,rasqual$pos,rasqual$pos,vcf_fn)
res = system(command,intern=TRUE)
header = str_split(res[1],'\t') 
header = str_extract(header[[1]],'(?<=])(.+?)(?=(:|$))')
dosage = str_split(res[2],'\t')[[1]]
dosage = data.table(sample = header,dosage)
dosage = dosage[sample!='ID',]
dosage[,dosage:=as.numeric(dosage)]

genotype_pc = read_genotype_pc(genotype_pc_fn)

merged = merge(dosage,genotype_pc,by='sample')
cor(merged[,2:5]) # Correlation between PC2 and dosage is 0.33826549
