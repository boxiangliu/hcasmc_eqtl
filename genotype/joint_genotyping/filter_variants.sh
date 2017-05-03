wd=../data/joint2
bcftools view -Oz -o $wd/recalibrated_variants.pass.vcf.gz -f PASS $wd/recalibrated_variants.vcf.gz
