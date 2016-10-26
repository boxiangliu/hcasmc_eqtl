in_vcf=$1
# in_vcf=../processed_data/160604_phasing/phased_and_imputed_gprobs/phased_and_imputed.chr6.vcf.gz

rename_vcf=${in_vcf/vcf.gz/rename.vcf.gz}
dr2_vcf=${rename_vcf/vcf.gz/dr2.vcf.gz}
hwe=${dr2_vcf/.vcf.gz/} # output file name will be appended .hwe
pass_hwe=${hwe}.pass_hwe.txt
hwe_vcf=${dr2_vcf/vcf.gz/hwe.vcf.gz}
indel_vcf=${hwe_vcf/vcf.gz/indellt51.vcf.gz}
rnasample_vcf=${indel_vcf/vcf.gz/rnasample.vcf.gz}
hg19_vcf=${rnasample_vcf/vcf.gz/hg19.vcf.gz}

rm $rename_vcf $dr2_vcf $hwe.hwe $hwe.log $pass_hwe $hwe_vcf $indel_vcf $indel_vcf.tbi $rnasample_vcf $hg19_vcf $hg19_vcf.tbi