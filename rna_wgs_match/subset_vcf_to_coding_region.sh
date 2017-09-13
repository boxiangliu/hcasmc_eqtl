out_dir=../processed_data/rna_wgs_match/subset_vcf_to_coding_region/
mkdir -p $out_dir
awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' ../data/gtex/gencode.v19.genes.v6p.hg19.bed | sed 's/chr//g' > $out_dir/coding_region.bed


bcftools view \
-R $out_dir/coding_region.bed \
-Oz -o $out_dir/chr22.vcf.gz \
../data/joint3/phased/phased_and_imputed.chr22.vcf.gz
