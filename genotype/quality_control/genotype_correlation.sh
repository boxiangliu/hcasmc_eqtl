mkdir -p ../processed_data/genotype/quality_control/genotype_correlation/

bcftools query -H \
-f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' \
../data/joint3/phased/phased_and_imputed.chr1.vcf.gz > \
../processed_data/genotype/quality_control/genotype_correlation/chr1_dosage.txt

Rscript genotype/quality_control/genotype_correlation.R