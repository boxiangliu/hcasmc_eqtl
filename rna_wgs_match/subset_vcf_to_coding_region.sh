out_dir=../processed_data/subset_vcf_to_coding_region/
mkdir -p $out_dir
awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' ../data/gtex/gencode.v19.genes.v6p.hg19.bed | sed 's/chr//g' > ../processed_data/subset_vcf_to_coding_region/coding_region.bed


bcftools view \
-R ../processed_data/subset_vcf_to_coding_region/coding_region.bed \
-Oz -o $out_dir/chr22.vcf.gz \
../data/joint3/phased/phased_and_imputed.chr22.vcf.gz
