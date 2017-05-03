# remove intermediate files:
wd=../data/joint2/
rm $wd/recalibrated_snps_raw_indels.vcf
rm $wd/recalibrated_snps_raw_indels.vcf.idx
mkdir $wd/variant_recalibration
mv $wd/recalibrate_SNP* $wd/variant_recalibration
mv $wd/recalibrate_INDEL* $wd/variant_recalibration
