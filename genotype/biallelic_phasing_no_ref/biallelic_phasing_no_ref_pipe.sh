# 16/06/28:
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts/gneotype/biallelic_phasing_no_ref/
processed_data=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/gneotype/biallelic_phasing_no_ref/
figures=/srv/persistent/bliu2/HCASMC_eQTL/figures/gneotype/biallelic_phasing_no_ref/
mkdir $scripts $processed_data $figures


# filter for biallelic variants:
wd=/srv/persistent/bliu2/HCASMC_eQTL/data/joint3
ln /srv/persistent/bliu2/HCASMC_eQTL/data/joint2/recalibrated_variants.pass.vcf.gz /srv/persistent/bliu2/HCASMC_eQTL/data/joint3/recalibrated_variants.pass.vcf.gz
bcftools view -m2 -M2 -Ov -o $wd/recalibrated_biallelic_variants.pass.vcf $wd/recalibrated_variants.pass.vcf.gz


# beagle imputation without reference: 
bash $scripts/run_beagle_impute.sh \
	/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/recalibrated_biallelic_variants.pass.vcf \
	/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/recalibrated_biallelic_variants.beagle


# post-imputation beagle QC:
cat /srv/persistent/bliu2/HCASMC_eQTL/data/joint3/recalibrated_biallelic_variants.beagle.vcf | \
	grep -v "^#" | \
	awk 'BEGIN {FS="\t|;|="; OFS="\t"; print "CHROM","POS","AR2","DR2","AF"} {print $1,$2,$9,$11,$13}' > \
	$processed_data/recalibrated_biallelic_SNP.r2.tsv
Rscript beagle_QC.R \
	-input=$processed_data/recalibrated_biallelic_SNP.r2.tsv \
	-figure_dir=$figures


# update sample names:
cd /srv/persistent/bliu2/HCASMC_eQTL/data/joint3
ln /srv/persistent/bliu2/HCASMC_eQTL/data/joint2/old_to_new_sample_name.txt /srv/persistent/bliu2/HCASMC_eQTL/data/joint3/old_to_new_sample_name.txt 
bcftools reheader -s old_to_new_sample_name.txt -o recalibrated_biallelic_variants.beagle.rename.vcf recalibrated_biallelic_variants.beagle.vcf


# filter for variants with dosage R2 >= 0.8:
cd /srv/persistent/bliu2/HCASMC_eQTL/data/joint3
bcftools view -e 'INFO/DR2<0.8' -o recalibrated_biallelic_variants.beagle.rename.dr2.vcf recalibrated_biallelic_variants.beagle.rename.vcf


# hwe filtering:
cd /srv/persistent/bliu2/HCASMC_eQTL/data/joint3
ln /srv/persistent/bliu2/HCASMC_eQTL/data/joint2/caucasian_for_hwe.txt /srv/persistent/bliu2/HCASMC_eQTL/data/joint3/caucasian_for_hwe.txt
grep -v "IMP" recalibrated_biallelic_variants.beagle.rename.dr2.vcf > recalibrated_biallelic_variants.beagle.rename.dr2.2.vcf
vcftools --vcf recalibrated_biallelic_variants.beagle.rename.dr2.2.vcf --keep caucasian_for_hwe.txt --hardy --out hwe_pval
rm recalibrated_biallelic_SNP.beagle.rename.dr2.2.vcf

# select sites with hwe > 1e-6:
tail -n +2 hwe_pval.hwe | awk 'BEGIN{OFS="\t"} {if ($6 >= 1e-6) print $1,$2}' > pass_hwe.txt
bcftools view -T pass_hwe.txt -Ov -o recalibrated_biallelic_variants.beagle.rename.dr2.hwe.vcf recalibrated_biallelic_variants.beagle.rename.dr2.vcf


# filter for MAF > 0.05:
cd /srv/persistent/bliu2/HCASMC_eQTL/data/joint3
bcftools view \
	--min-af 0.05 --max-af 0.95 \
	-Ov -o recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf \
	recalibrated_biallelic_variants.beagle.rename.dr2.hwe.vcf


# change ID: 
cd /srv/persistent/bliu2/HCASMC_eQTL/data/joint3
bcftools annotate -Oz \
	-o recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf.id.gz \
	--set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
	recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf
tabix -p vcf recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf.id.gz

