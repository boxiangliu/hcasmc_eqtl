wd=../data/joint2/

# change chromosome names from "chr*" to "*":
for i in $(seq 1 22) X Y M; do 
if [[ $i=="M" ]]; then 
	echo "chrM MT" >> $wd/hg19_to_GRCh37.txt
else 
	echo "chr$i $i" >> $wd/hg19_to_GRCh37.txt
fi
done
bcftools annotate --rename-chrs $wd/hg19_to_GRCh37.txt -Oz -o $wd/recalibrated_biallelic_SNP.pass.GRCh37.vcf.gz $wd/recalibrated_biallelic_SNP.pass.vcf.gz
