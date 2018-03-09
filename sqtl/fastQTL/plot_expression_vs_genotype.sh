out_dir=../processed_data/sqtl/fastQTL/plot_expression_vs_genotype/
mkdir -p $out_dir

# Convert VCF file to dosage:
bcftools query \
-f '%ID[\t%DS]\n' \
-r 19:44153100-44153100 \
-H \
../data/joint3/phased/phased_and_imputed.chr19.vcf.gz | \
sed 's/\[[0-9]\+\]//g' | \
sed 's/# //g' | \
sed 's/:DS//g' \
> $out_dir/dosage.tsv


# Subset expression file:
grep -e 19:44244370:44248924:clu_12328 -e Chr \
../data/rnaseq2/leafcutter_wasp/cluster/sqtl_perind.counts.gz.phen_chr19 | \
cut -f4- \
> $out_dir/expression.tsv


# Make plot:
Rscript sqtl/fastQTL/plot_expression_vs_genotype.R 