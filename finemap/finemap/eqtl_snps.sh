concat_vcf=../processed_data/finemap/finemap/rasqual_data/all.vcf.gz
bcftools concat ../data/joint3/asvcf/phased_and_imputed.chr{1..22}.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz -Oz -o $concat_vcf
tabix -p vcf $concat_vcf
eqtl_vcf=../processed_data/finemap/finemap/rasqual_data/eqtl.vcf
bcftools view -R ../processed_data/finemap/finemap/rasqual_data/eqtl.chrpos.txt $concat_vcf > $eqtl_vcf


concat_vcf=../processed_data/finemap/finemap/rasqual_data/all.id.vcf.gz
bcftools concat ../data/joint3/asvcf_sid/phased_and_imputed.chr{1..22}.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz -Oz -o $concat_vcf
tabix -p vcf $concat_vcf
eqtl_vcf=../processed_data/finemap/finemap/rasqual_data/eqtl.id.vcf
bcftools view -R ../processed_data/finemap/finemap/rasqual_data/eqtl.chrpos.txt $concat_vcf | vcf-sort -c | uniq > $eqtl_vcf
