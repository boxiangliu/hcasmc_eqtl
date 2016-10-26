in_vcf=$1
# in_vcf=../processed_data/160604_phasing/phased_and_imputed_gprobs/phased_and_imputed.chr6.vcf.gz


# variables: 
tools=/srv/persistent/bliu2/tools
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts/
rename_vcf=${in_vcf/vcf.gz/rename.vcf.gz}
dr2_vcf=${rename_vcf/vcf.gz/dr2.vcf.gz}
hwe=${dr2_vcf/.vcf.gz/} # output file name will be appended .hwe
pass_hwe=${hwe}.pass_hwe.txt
hwe_vcf=${dr2_vcf/vcf.gz/hwe.vcf.gz}
indel_vcf=${hwe_vcf/vcf.gz/indellt51.vcf.gz}
rnasample_vcf=${indel_vcf/vcf.gz/rnasample.vcf.gz}
hg19_vcf=${rnasample_vcf/vcf.gz/hg19.vcf.gz}


# rename the sample columns:
echo "renaming sampling columns..."
bcftools reheader -s ../data/joint3/old_to_new_sample_name.txt -o $rename_vcf $in_vcf


# select DR2 greater than 0.8: 
echo "filtering DR2<0.8..."
bcftools view -e 'INFO/DR2<0.8' -Oz -o $dr2_vcf $rename_vcf


# select sites with hwe > 1e-6:
echo "filtering HWE<1e-6..."
vcftools --gzvcf $dr2_vcf --keep ../data/joint3/caucasian_for_hwe.txt --hardy --out $hwe
tail -n +2 $hwe.hwe | awk 'BEGIN{OFS="\t"} {if ($6 >= 1e-6) print $1,$2}' > ${pass_hwe}
bcftools view -T ${pass_hwe} -Oz -o ${hwe_vcf} ${dr2_vcf}


# filter out INDEL greater than 51bp:
echo "filtering indel < 51bp..."
zcat ${hwe_vcf} | awk '{if (length($4)<=51 && length($5)<=51) print $0}' | bgzip > ${indel_vcf}


# only keep VCF sample with RNAseq data:
echo "filtering RNAseq samples in VCF file..."
tabix -p vcf ${indel_vcf}
bcftools view \
	-S ../processed_data/160604_phasing/phased_and_imputed_gprobs/sample_list.txt \
	-Oz -o ${rnasample_vcf} ${indel_vcf}


# convert VCF from GRCh37 to hg19:
echo "converting GRCh37 to hg19..."
bcftools annotate --rename-chrs ../data/joint2/GRCh37_to_hg19.txt \
	-Oz -o ${hg19_vcf} ${rnasample_vcf}

# index:
echo "indexing VCF..."
tabix -p vcf ${hg19_vcf}


# make ASVCF: 
# ls ../data/rnaseq2/wasp/*/wasp.keep.merged.rmdup.sorted.bam > ../processed_data/rasqual/bam.list.txt
echo "counting ASE..."
bash /srv/persistent/bliu2/tools/rasqual/src/ASVCF/createASVCF.sh \
	../processed_data/rasqual/bam.list.txt \
	${hg19_vcf}