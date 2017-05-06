scripts=genotype/phasing_with_1kg/
processed_data=../processed_data/genotype/phasing_with_1kg/

# convert hg19 coordinate to GRCh37 coordinate:
bash $scripts/convert_chrom_names.sh ../data/joint2/recalibrated_biallelic_SNP.beagle.rename.dr2.hwe.vcf.gz ../data/joint2/recalibrated_biallelic_SNP.beagle.rename.dr2.hwe.GRCh37.vcf.gz ../data/joint2/hg19_to_GRCh37.txt
bash $scripts/convert_chrom_names.sh ../data/joint3/recalibrated_variants.pass.vcf.gz ../data/joint3/recalibrated_variants.pass.GRCh37.vcf.gz ../data/joint2/hg19_to_GRCh37.txt


# make a list of non-european samples:
Rscript $scripts/make_list_of_non_caucasian_samples.hcasmc.R ../processed_data/rna_wgs_match.reduced_050616.xlsx $processed_data/non_caucasian.hcasmc.txt
bash $scripts/make_list_of_non_EUR_samples.1kg.sh /srv/persistent/bliu2/tools/beagle/reference/integrated_call_samples_v3.20130502.ALL.panel $processed_data/non_caucasian.1kg.txt
cat $processed_data/non_caucasian.1kg.txt $processed_data/non_caucasian.hcasmc.txt > $processed_data/non_caucasian.txt


# run comform-gt on filtered variants:
for i in {1..22} X; do
	bash $scripts/conform_gt.core.sh \
	/srv/persistent/bliu2/tools/beagle/reference/chr$i.1kg.phase3.v5a.vcf.gz \
	../data/joint2/recalibrated_biallelic_SNP.beagle.rename.dr2.hwe.GRCh37.vcf.gz \
	$i \
	$processed_data/conform_gt/mod.chr$i \
	$processed_data/non_caucasian.txt &
done


# run comform-gt on all recalibrated variants: 
parallel -j6 bash $scripts/conform_gt.core.sh \
	/srv/persistent/bliu2/tools/beagle/reference/chr{}.1kg.phase3.v5a.vcf.gz \
	/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/recalibrated_variants.pass.GRCh37.vcf.gz \
	{} \
	$processed_data/conform_gt_all_variants/mod.chr{} \
	$processed_data/non_caucasian.txt ::: {1..22} X
