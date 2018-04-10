library(data.table)
library(stringr)

DNA_vcf_fn = '../data/joint3/orig/recalibrated_variants.pass.vcf.gz'
RNA_vcf_fn = '../data/joint3/asvcf/phased_and_imputed.chrX.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz'
old_to_new_fn = '../data/joint3/orig/old_to_new_sample_name.txt'

get_VCF_samples = function(vcf_fn){
	command = sprintf('bcftools view -h %s | tail -n1',vcf_fn)
	print(command)
	result = system(command,intern=TRUE)
	split_result = str_split(result,'\t')[[1]]
	samples = split_result[10:length(split_result)]
	return(samples)
}

update_sample_names = function(sample_names,old_to_new_fn){
	old_to_new = fread(old_to_new_fn)
	setnames(old_to_new,c('old','new'))
	sample_names[sample_names %in% old_to_new$old] = old_to_new$new
	return(sample_names)
}

DNA_samples = get_VCF_samples(DNA_vcf_fn)
DNA_samples = update_sample_names(DNA_samples,old_to_new_fn)
RNA_samples = get_VCF_samples(RNA_vcf_fn)

DNA_samples[!(DNA_samples %in% RNA_samples)]
# "1497"  "1923"  "2115"  "2161"  "24156" "2477"