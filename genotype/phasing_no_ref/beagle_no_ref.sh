wd=../data/joint2/
java=/srv/persistent/bliu2/tools/jre1.8.0_91/bin/java
beagle=/srv/persistent/bliu2/tools/beagle/beagle.03May16.862.jar

# impute and phase:
for i in $(seq 1 22);do
	$java -Xmx4096m -jar $beagle nthreads=2 chrom=chr$i gl=$wd/recalibrated_biallelic_SNP.pass.vcf out=$wd/recalibrated_biallelic_SNP.beagle.chr$i &
	tabix -p vcf $wd/recalibrated_biallelic_SNP.beagle.chr$i.vcf.gz
done


# concatenate all chromosomes:
bcftools concat -Oz -o recalibrated_biallelic_SNP.beagle.vcf.gz recalibrated_biallelic_SNP.beagle.chr{1..22}.vcf.gz
# I made sure the beagle output is complete by comparing the CHROM and POS fields of 
# recalibrated_biallelic_SNP.beagle.vcf and recalibrated_biallelic_SNP.pass.vcf

# move beagle intermediate files to folder: 
mkdir $wd/beagle_no_ref
mv $wd/recalibrated_biallelic_SNP.beagle.chr* $wd/beagle_no_ref

