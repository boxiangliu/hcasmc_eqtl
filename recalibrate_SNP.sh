#!/bin/bash 
# bosh liu
# durga
# 2016/05/12
# run GATK variantRecalibrator

# command args: 
input=$1


# paths: 
GATK=/usr/bin/GenomeAnalysisTK.jar
HG19=/srv/persistent/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
GATK_BUNDLE=/srv/persistent/bliu2/shared/gatk_bundle_2.8_hg19

wd=../data/joint2


# change working directory: 
cd $wd 


# run VariantRecalibrator:
# follows https://www.broadinstitute.org/gatk/guide/article?id=1259
java -Xmx16g -jar $GATK -nt 4 -T VariantRecalibrator -R $HG19 -input $input \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATK_BUNDLE/hapmap_3.3.hg19.sites.vcf \
	-resource:omni,known=false,training=true,truth=true,prior=12.0 $GATK_BUNDLE/1000G_omni2.5.hg19.sites.vcf \
	-resource:1000G,known=false,training=true,truth=false,prior=10.0 $GATK_BUNDLE/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATK_BUNDLE/dbsnp_138.hg19.vcf \
	-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff \
	-mode SNP \
	-tranche 100.0 -tranche 99.9 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.1 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 90.0 \
	-recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches \
	-rscriptFile recalibrate_SNP_plots.R


# run ApplyRecalibration:
java -Xmx3g -jar $GATK -T ApplyRecalibration -R $HG19 -input $input \
	-mode SNP \
	--ts_filter_level 99.5 \
	-recalFile recalibrate_SNP.recal \
	-tranchesFile recalibrate_SNP.tranches \
	-o recalibrated_snps_raw_indels.vcf