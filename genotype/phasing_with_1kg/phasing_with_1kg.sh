#### phasing 
# 16/06/04:
# setup:
scripts=genotype/phasing_with_1kg
processed_data=../processed_data/genotype/phasing_with_1kg/
mkdir $scripts $processed_data 



# Beagle phasing with reference no imputation:
mkdir $processed_data/phased/
n=0
for i in $(seq 1 22);do
bash $scripts/run_beagle_phasing.sh \
	$processed_data/conform_gt/mod.chr$i.vcf.gz \
	$processed_data/phased/phased.chr$i \
	/srv/persistent/bliu2/tools/beagle/reference/chr$i.1kg.phase3.v5a.vcf.gz \
	/srv/persistent/bliu2/tools/beagle/reference/plink.chr$i.GRCh37.map \
	$i &
n=$(($n+1))
if [[ $n -gt 10 ]]; then wait; n=0; fi
done

# merge vcf: 
bcftools concat -Oz -o $processed_data/phased/phased.vcf.gz $processed_data/phased/phased.{1..22}.vcf.gz



# Beagle phasing with reference with imputation (regular output):
mkdir ../processed_data/160604_phasing/phased_and_imputed
parallel -j6 /srv/persistent/bliu2/tools/jre1.8.0_91/bin/java -Xmx8g -jar \
	/srv/persistent/bliu2/tools/beagle/beagle.27Jul16.86a.jar \
	nthreads=4 \
	chrom={} \
	gt=../processed_data/160604_phasing/conform_gt_all_variants/mod.chr{}.vcf.gz \
	out=../processed_data/160604_phasing/phased_and_imputed/phased_and_imputed.chr{} \
	ref=/srv/persistent/bliu2/tools/beagle/reference/chr{}.1kg.phase3.v5a.bref \
	map=/srv/persistent/bliu2/tools/beagle/reference/plink.chr{}.GRCh37.map \
	impute=true ::: {1..22} X


# Beagle phasing with reference with imputation (output genotype probability):
mkdir ../processed_data/160604_phasing/phased_and_imputed_gprobs
parallel -j12 /srv/persistent/bliu2/tools/jre1.8.0_91/bin/java -Xmx8g -jar \
	/srv/persistent/bliu2/tools/beagle/beagle.27Jul16.86a.jar \
	nthreads=2 \
	chrom={} \
	gt=../processed_data/160604_phasing/conform_gt_all_variants/mod.chr{}.vcf.gz \
	out=../processed_data/160604_phasing/phased_and_imputed_gprobs/phased_and_imputed.chr{} \
	ref=/srv/persistent/bliu2/tools/beagle/reference/chr{}.1kg.phase3.v5a.bref \
	map=/srv/persistent/bliu2/tools/beagle/reference/plink.chr{}.GRCh37.map \
	impute=true gprobs=true ::: {1..22} X

# link result to data directory: 
ln -s /srv/persistent/bliu2/HCASMC_eQTL/processed_data/genotype/phasing_with_1kg/phased_and_imputed_gprobs /srv/persistent/bliu2/HCASMC_eQTL/data/joint3/phased_imputed


# Beagle phasing without reference: 
mkdir $processed_data/phased_no_ref/
n=0
for i in $(seq 1 22);do
bash $scripts/run_beagle_phasing_without_reference.sh \
	../data/joint2/recalibrated_biallelic_SNP.beagle.rename.dr2.hwe.GRCh37.vcf.gz \
	$processed_data/phased_no_ref/phased_no_ref.chr$i \
	/srv/persistent/bliu2/tools/beagle/reference/plink.chr$i.GRCh37.map \
	$i &
n=$(($n+1))
if [[ $n -gt 10 ]]; then wait; n=0; fi
done

# merge vcf:
for i in {1..22}; do tabix -p vcf $processed_data/phased_no_ref/phased_no_ref.chr$i.vcf.gz; done
bcftools concat -Oz -o $processed_data/phased_no_ref/phased_no_ref.vcf.gz $processed_data/phased_no_ref/phased_no_ref.chr{1..22}.vcf.gz

# convert GRCh37 to hg19: 
bash $scripts/convert_chrom_names.sh \
	/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160604_phasing/phased_no_ref/phased_no_ref.vcf.gz \
	/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160604_phasing/phased_no_ref/phased_no_ref.hg19.vcf \
	/srv/persistent/bliu2/HCASMC_eQTL/data/joint2/GRCh37_to_hg19.txt
