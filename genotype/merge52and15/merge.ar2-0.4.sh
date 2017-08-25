wd=/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/
mkdir -p $wd/67total_ar2-0.4

for i in {1..22} X; do

echo $i

bcftools reheader \
-s /srv/persistent/bliu2/HCASMC_eQTL/data/joint2/old_to_new_sample_name.txt \
$wd/phased/phased_and_imputed.chr$i.vcf.gz -o $wd/67total_ar2-0.4/chr$i.52.temp.vcf.gz

tabix -p vcf $wd/67total_ar2-0.4/chr$i.52.temp.vcf.gz

bcftools view \
-i "INFO/AR2>=0.4" \
-S /srv/persistent/bliu2/HCASMC_eQTL/processed_data/genotype/phasing_with_1kg/phased_and_imputed_gprobs/sample_list.txt \
-Oz -o $wd/67total_ar2-0.4/chr$i.52.vcf.gz $wd/67total_ar2-0.4/chr$i.52.temp.vcf.gz

tabix -p vcf $wd/67total_ar2-0.4/chr$i.52.vcf.gz


ln $wd/15new/beagle/chr$i.ar2-0.4.vcf.gz \
$wd/67total_ar2-0.4/chr$i.15.vcf.gz

tabix -p vcf $wd/67total_ar2-0.4/chr$i.15.vcf.gz

bcftools merge -i AR2:join,DR2:join \
$wd/67total_ar2-0.4/chr$i.52.vcf.gz \
$wd/67total_ar2-0.4/chr$i.15.vcf.gz -Oz -o $wd/67total_ar2-0.4/chr$i.vcf.gz

tabix -p vcf $wd/67total_ar2-0.4/chr$i.vcf.gz

bcftools view -g ^miss \
$wd/67total_ar2-0.4/chr$i.vcf.gz \
-Oz -o $wd/67total_ar2-0.4/chr$i.vcf.comm.gz 

tabix -p vcf $wd/67total_ar2-0.4/chr$i.vcf.comm.gz

done