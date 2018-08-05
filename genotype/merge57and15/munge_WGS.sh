out_dir=/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/72total/57WGS_ar2-0.4/

for i in {1..22} X; do
bcftools view \
--include "INFO/AR2 >= 0.4" \
--samples ^24156 \
/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/phased/phased_and_imputed.chr$i.vcf.gz \
-Oz -o $out_dir/chr$i.toreheader.vcf.gz

bcftools reheader \
--samples /srv/persistent/bliu2/HCASMC_eQTL/scripts/genotype/merge57and15/WGS_sample_names.txt \
$out_dir/chr$i.toreheader.vcf.gz \
-o $out_dir/chr$i.vcf.gz

tabix $out_dir/chr$i.vcf.gz
rm $out_dir/chr$i.toreheader.vcf.gz
done