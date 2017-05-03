#!/bin/bash
# bosh liu
# durga
# detect sample contamination using model-based method: 

# command arg: 
sample=$1
ethnicity=$2 # should be one of EUR, AFR, EAS, AMR or SAS.


# paths: 
processed_data=../processed_data/genotype/quality_control/detect_WGS_contamination/


# overwrite 1000G VCF AF field with ${ethnicity}_AF:
if [[ ! -f $processed_data/chr20.${ethnicity}_AF.vcf.gz ]];then
	bcftools view -Ou -s . --force-samples --no-update --types snps /srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | bcftools annotate -Ov -x ^INFO/${ethnicity}_AF | sed -e "s/${ethnicity}_AF/AF/" | bcftools view -Oz -o $processed_data/chr20.${ethnicity}_AF.vcf.gz
fi

# subset to chr20 and 
# modify chromosome labels from hg19 format to GRCh37 format: 
samtools view -h /mnt/data/WGS_HCASMC/$sample/recal_reads.bam chr20 | sed 's/chr20/20/' | samtools view -h -b -o $processed_data/recal_reads.$sample.chr20.bam -
samtools index $processed_data/recal_reads.$sample.chr20.bam

# run verfyBAMID: 
verifyBamID --vcf $processed_data/chr20.${ethnicity}_AF.vcf.gz --bam $processed_data/recal_reads.$sample.chr20.bam --chip-none --precise --verbose --noPhoneHome --out $processed_data/verifyBAMID.$sample



