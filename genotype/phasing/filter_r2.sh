cd ../data/joint2/

# filter for variants with dosage R2 >= 0.8:
bcftools view -e 'INFO/DR2<0.8' -o recalibrated_biallelic_SNP.beagle.rename.dr2.vcf.gz recalibrated_biallelic_SNP.beagle.rename.vcf.gz

