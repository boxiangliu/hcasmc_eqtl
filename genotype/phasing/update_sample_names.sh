# update sample names:
wd=../data/joint2/
cd $wd
echo "CA1401 1401" >> old_to_new_sample_name.txt
echo "2102 2105" >> old_to_new_sample_name.txt
echo "2109 1508" >> old_to_new_sample_name.txt
echo "289727 2999" >> old_to_new_sample_name.txt
echo "313605 317155" >> old_to_new_sample_name.txt
bcftools reheader -s old_to_new_sample_name.txt -o recalibrated_biallelic_SNP.beagle.rename.vcf.gz recalibrated_biallelic_SNP.beagle.vcf.gz
