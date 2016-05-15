#!/bin/bash 
# bosh liu
# durga
# 2016/05/12
# run genotypeGVCF

# command args:
chrom=$1


# paths: 
GATK=/usr/bin/GenomeAnalysisTK.jar
HG19=/srv/persistent/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
GATK_BUNDLE=/srv/persistent/bliu2/shared/gatk_bundle_2.8_hg19

wd=../data/joint2
sample_list=sample_list.txt
out_vcf=raw_variants.$chrom.vcf


# change working directory: 
cd $wd 


# build commands:
cmd="java -Xmx4g -jar $GATK -T GenotypeGVCFs -R $HG19 --dbsnp $GATK_BUNDLE/dbsnp_138.hg19.vcf -L $chrom"

while read sample; do
cmd="$cmd --variant /mnt/data/WGS_HCASMC/$sample/raw_variants.g.vcf.gz"
done < $sample_list

cmd="$cmd -o $out_vcf"


# run command:
eval $cmd