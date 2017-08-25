wd=/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/

for i in {1..22} X; do

echo $i

bcftools annotate \
-x INFO/AF,INFO/AC,INFO/AN,FORMAT/AS \
$wd/asvcf/phased_and_imputed.chr$i.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz \
-Oz -o $wd/67total/chr$i.52.vcf.gz

tabix -p vcf $wd/67total/chr$i.52.vcf.gz


bcftools annotate \
--rename-chrs /srv/persistent/bliu2/HCASMC_eQTL/scripts/genotype/merge52and15/old2new_chrom_name.txt \
$wd/15new/beagle/chr$i.ar2-0.4.vcf.gz \
-Oz -o $wd/67total/chr$i.15.vcf.gz

tabix -p vcf $wd/67total/chr$i.15.vcf.gz

bcftools merge -i AR2:join,DR2:join \
$wd/67total/chr$i.52.vcf.gz \
$wd/67total/chr$i.15.vcf.gz -Oz -o $wd/67total/chr$i.vcf.gz

tabix -p vcf $wd/67total/chr$i.vcf.gz

bcftools view -g ^miss \
$wd/67total/chr$i.vcf.gz \
-Oz -o $wd/67total/chr$i.vcf.comm.gz 

tabix -p vcf $wd/67total/chr$i.vcf.comm.gz

done