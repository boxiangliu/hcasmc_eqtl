wgs_dir=/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/72total/57WGS_ar2-0.4/
array_dir=/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/15new/beagle/

out_dir=/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/72total/merged/
mkdir -p $out_dir

for i in {1..22} X; do
bcftools annotate \
-x INFO/AF,INFO/AC,INFO/AN \
$wgs_dir/chr$i.vcf.gz \
-Oz -o $out_dir/chr$i.57wgs.vcf.gz &
done
wait 

for i in {1..22} X; do
tabix -p vcf $out_dir/chr$i.57wgs.vcf.gz &
done 
wait 

for i in {1..22} X; do
bcftools merge -i AR2:join,DR2:join \
$out_dir/chr$i.57wgs.vcf.gz \
$array_dir/chr$i.ar2-0.4.vcf.gz \
-Oz -o $out_dir/chr$i.vcf.gz &
done
wait

