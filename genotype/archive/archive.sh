# joint genotype calling: 
/mnt/data/WGS_HCASMC/joint/recalibrated_variants.vcf:
	python genotype_gvcfs.py /mnt/data/WGS_HCASMC/joint/ ../sample_list.txt 2> log/genotype_gvcfs.20160319.log 
	# relocate temporary files:
	cd /mnt/data/WGS_HCASMC/joint
	mkdir genotype_gvcfs_temp
	mv recalibrated_snps_raw_indels.vcf* genotype_gvcfs_temp
	mv recalibrated_variants.vcf* genotype_gvcfs_temp
	mv recalibrate_INDEL* genotype_gvcfs_temp
	mv recalibrate_SNP* genotype_gvcfs_temp


# genotype refinement: 
/mnt/data/WGS_HCASMC/joint/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf.gz: /mnt/data/WGS_HCASMC/joint/recalibrated_variants.vcf.gz
	python genotype_refinement.py /mnt/data/WGS_HCASMC/joint/ recalibrated_variants.GRCh37.vcf.gz 2>> log/genotype_refinement.log
	cd /mnt/data/WGS_HCASMC/joint/; mkdir genotype_refinement_temp
	mv recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf.idx genotype_refinement_temp
	mv recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf genotype_refinement_temp
	mv recalibrated_variants.GRCh37.postCGP.vcf.idx genotype_refinement_temp
	mv recalibrated_variants.GRCh37.postCGP.vcf genotype_refinement_temp


# vcf quality control using QC3: 
input4 = /mnt/data/WGS_HCASMC/joint/genotype_refinement_temp/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf
output4 = /mnt/data/WGS_HCASMC/joint/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf.QC3
$(output4): $(input4)
	perl /srv/persistent/bliu2/tools/QC3/qc3.pl -m v -i $(input4) -a /srv/persistent/bliu2/tools/annovar/humandb -o $(output4)


input5 = /mnt/data/WGS_HCASMC/joint/genotype_gvcfs_temp/recalibrated_variants.vcf
output5 = /mnt/data/WGS_HCASMC/joint/recalibrated_variants.vcf.QC3
$(output5): $(input5)
	perl /srv/persistent/bliu2/tools/QC3/qc3.pl -m v -i $(input5) -a /srv/persistent/bliu2/tools/annovar/humandb -o $(output5)


# vcf quality control using plink/seq: 
../processed_data/recalibrated_variants.vcf.pseq.vstats: /mnt/data/WGS_HCASMC/joint/genotype_gvcfs_temp/recalibrated_variants.vcf
	pseq /mnt/data/WGS_HCASMC/joint/genotype_gvcfs_temp/recalibrated_variants.vcf v-stats --mask reg=chr22 any.filter.ex > ../processed_data/recalibrated_variants.vcf.pseq.vstats
	pseq /mnt/data/WGS_HCASMC/joint/genotype_gvcfs_temp/recalibrated_variants.vcf i-stats --mask reg=chr22 any.filter.ex > ../processed_data/recalibrated_variants.vcf.pseq.istats
	pseq /mnt/data/WGS_HCASMC/joint/genotype_refinement_temp/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf v-stats --mask reg=chr22 any.filter.ex > ../processed_data/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf.pseq.vstats
	pseq /mnt/data/WGS_HCASMC/joint/genotype_refinement_temp/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf i-stats --mask reg=chr22 any.filter.ex > ../processed_data/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf.pseq.istats



# plot plink/seq istats files:
input6 = ../processed_data/recalibrated_variants.vcf.pseq.istats
figure6 = ../figures/pseq_istats_recalibrated_variants.pdf
$(figure6): $(input6)
	Rscript plot_pseq_stats.R ../processed_data

# remove multiallelic loci: 
input20 = ../data/joint/recalibrated_variants.GRCh37.vcf.gz
output20 = ../data/joint/recalibrated_variants.GRCh37.biallelic.vcf.gz
$(output20): $(input20)
	# Keeping only biallelic SNPs and Indels 
	bcftools view -m2 -M2 -v snps,indels $(input20) -Oz -o $(output20)
	tabix -p vcf $(output20)

	# what percentage are multi-allelic sites: 
	bcftools view -H $(input20) | wc -l  # 26604334
	bcftools view -H $(output20) | wc -l # 25275762
	# a=26604334
	# b=25275762
	# (a-b)/a = 0.04993818

# remove multiallelic loci and missing genotypes: 
input10 = ../data/joint/recalibrated_variants.GRCh37.vcf.gz
output10 = $(input10:vcf.gz=biallelic.nomissing.vcf.gz)
$(output10): $(input10)
	# Keeping only biallelic SNPs and Indels; removing missing genotypes. 
	bcftools view -m2 -M2 -v snps,indels -g ^miss $(input10) -Oz -o $(output10)
	tabix -p vcf $(output10)


