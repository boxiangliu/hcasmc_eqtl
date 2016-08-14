import vcf,sys
vcf_file='/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf'
vcf_reader = vcf.Reader(open(vcf_file,'r'))
for sample in vcf_reader.samples:
	print(sample)