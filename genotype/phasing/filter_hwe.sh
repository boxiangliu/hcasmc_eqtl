# subset for Caucasian individuals to apply HWE filter:

Rscript genotype/phasing/caucasian_individual_for_hwe.R

# hwe filtering:
cd ../data/joint2/
zgrep -v "IMP" recalibrated_biallelic_SNP.beagle.rename.dr2.vcf.gz > recalibrated_biallelic_SNP.beagle.rename.dr2.2.vcf
vcftools --vcf recalibrated_biallelic_SNP.beagle.rename.dr2.2.vcf --keep caucasian_for_hwe.txt --hardy --out hwe_pval
rm recalibrated_biallelic_SNP.beagle.rename.dr2.2.vcf

# select sites with hwe > 1e-6:
tail -n +2 hwe_pval.hwe | awk 'BEGIN{OFS="\t"} {if ($6 >= 1e-6) print $1,$2}' > pass_hwe.txt
bcftools view -T pass_hwe.txt -Ov -o recalibrated_biallelic_SNP.beagle.rename.dr2.hwe.vcf.gz recalibrated_biallelic_SNP.beagle.rename.dr2.vcf.gz
