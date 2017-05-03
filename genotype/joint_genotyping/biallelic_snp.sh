wc=../data/joint2/
bcftools view -m2 -M2 -v snps -Oz -o $wd/recalibrated_biallelic_SNP.pass.vcf.gz $wd/recalibrated_variants.pass.vcf.gz 
