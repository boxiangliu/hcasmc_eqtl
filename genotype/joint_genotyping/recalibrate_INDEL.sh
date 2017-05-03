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
	-resource:mills,known=false,training=true,truth=true,prior=12.0 $GATK_BUNDLE/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATK_BUNDLE/dbsnp_138.hg19.vcf \
	-an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff \
	-mode INDEL \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
	--maxGaussians 4 \
	-recalFile recalibrate_INDEL.recal \
	-tranchesFile recalibrate_INDEL.tranches \
	-rscriptFile recalibrate_INDEL_plots.R



# run ApplyRecalibration:
java -Xmx4g -jar $GATK -T ApplyRecalibration -R $HG19 -input $input \
	-mode INDEL \
	--ts_filter_level 98.0 \
	-recalFile recalibrate_INDEL.recal \
	-tranchesFile recalibrate_INDEL.tranches \
	-o recalibrated_variants.vcf