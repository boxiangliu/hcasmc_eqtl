# extract dosage field:
# the output column order will be the same as that in sample_list.txt
wd=../data/joint2/
bcftools query -S $wd/sample_list.txt -f '%CHROM\_%POS\_%REF\_%ALT[\t%DS]\n' -o ../processed_data/genotype/phasing/extract_dosage/dosage.tsv $wd/recalibrated_biallelic_SNP.beagle.rename.vcf.gz
bcftools query -S $wd/sample_list.txt -f '%CHROM\_%POS\_%REF\_%ALT[\t%GT]\n' -o ../processed_data/genotype/phasing/genotype.tsv $wd/recalibrated_biallelic_SNP.beagle.rename.vcf.gz
bcftools query -S $wd/sample_list.txt -f '%CHROM\_%POS\_%REF\_%ALT[\t%GP]\n' -o ../processed_data/genotype/phasing/genotype_probability.tsv $wd/recalibrated_biallelic_SNP.beagle.rename.vcf.gz
