bcftools reheader -s ../data/joint3/old_to_new_sample_name.txt \
	-o ../processed_data/160604_phasing/phased_and_imputed/phased_and_imputed.chr6.rename.vcf.gz \
	../processed_data/160604_phasing/phased_and_imputed/phased_and_imputed.chr6.vcf.gz 

# select DR2 greater than 0.8: 
bcftools view -e 'INFO/DR2<0.8' -Oz -o ../processed_data/160604_phasing/phased_and_imputed/phased_and_imputed.chr6.rename.dr2.vcf.gz ../processed_data/160604_phasing/phased_and_imputed/phased_and_imputed.chr6.rename.vcf.gz


# select sites with hwe > 1e-6:
vcftools --gzvcf ../processed_data/160604_phasing/phased_and_imputed/phased_and_imputed.chr6.rename.dr2.vcf.gz \
	--keep ../data/joint3/caucasian_for_hwe.txt --hardy \
	--out ../processed_data/160604_phasing/phased_and_imputed/hwe_pval
tail -n +2 ../processed_data/160604_phasing/phased_and_imputed/hwe_pval.hwe | awk 'BEGIN{OFS="\t"} {if ($6 >= 1e-6) print $1,$2}' > ../processed_data/160604_phasing/phased_and_imputed/pass_hwe.txt
bcftools view -T ../processed_data/160604_phasing/phased_and_imputed/pass_hwe.txt \
	-Oz -o ../processed_data/160604_phasing/phased_and_imputed/phased_and_imputed.chr6.rename.dr2.hwe.vcf.gz \
	../processed_data/160604_phasing/phased_and_imputed/phased_and_imputed.chr6.rename.dr2.vcf.gz


# filter out INDEL greater than 51bp:
zcat ../processed_data/160604_phasing/phased_and_imputed/phased_and_imputed.chr6.rename.dr2.vcf.gz | \
	awk '{if (length($4)<=51 && length($5)<=51) print $0}' | \
	bgzip > ../processed_data/160604_phasing/phased_and_imputed/phased_and_imputed.chr6.rename.dr2.indellt51.vcf.gz


# 
cd /srv/persistent/bliu2/HCASMC_eQTL/data//rnaseq2/wasp/
ls -d */ | sed -e "s/\///" > $processed_data/160604_phasing/phased_and_imputed/sample_list.txt
cd $scripts
tabix -p vcf ../processed_data/160604_phasing/phased_and_imputed/phased_and_imputed.chr6.rename.dr2.indellt51.vcf.gz
bcftools view -S ../processed_data/160604_phasing/phased_and_imputed/sample_list.txt -Oz -o ../processed_data/160604_phasing/phased_and_imputed/phased_and_imputed.chr6.rename.dr2.indellt51.rnasamples.vcf.gz ../processed_data/160604_phasing/phased_and_imputed/phased_and_imputed.chr6.rename.dr2.indellt51.vcf.gz
ls /srv/persistent/bliu2/HCASMC_eQTL/data//rnaseq2/wasp/*/wasp.keep.merged.rmdup.sorted.bam > ../processed_data/rasqual/bam.list.txt

bcftools annotate --rename-chrs /srv/persistent/bliu2/HCASMC_eQTL/data/joint2/GRCh37_to_hg19.txt \
	-Oz -o ../processed_data/160604_phasing/phased_and_imputed/phased_and_imputed.chr6.rename.dr2.indellt51.rnasamples.hg19.vcf.gz \
	../processed_data/160604_phasing/phased_and_imputed/phased_and_imputed.chr6.rename.dr2.indellt51.rnasamples.vcf.gz


tabix -p vcf ../processed_data/160604_phasing/phased_and_imputed/phased_and_imputed.chr6.rename.dr2.indellt51.rnasamples.hg19.vcf.gz
# make ASVCF: 
cd $tools/rasqual
bash /srv/persistent/bliu2/tools/rasqual/src/ASVCF/createASVCF.sh $processed_data/rasqual/bam.list.txt $processed_data/160604_phasing/phased_and_imputed/phased_and_imputed.chr6.rename.dr2.indellt51.rnasamples.hg19.vcf.gz