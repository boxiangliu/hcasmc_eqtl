#!/bin/bash 
# boxiang liu
# durga
# calculate LD

indir=$1
outdir=$2
tkg=/srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.concat.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf
target_variants=$outdir/target_variants.txt
target_vcf=$outdir/target_variants.vcf

# make variant list:
cat $indir/*.txt | awk 'BEGIN{OFS="\t"}{print $1,$2}' | uniq > $target_variants


# subset to listed variants:
bcftools view \
-R $target_variants \
-Ov -o $target_vcf \
/srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.concat.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz


# generate sorted bed: 
plink \
--vcf $target_vcf --keep-allele-order \
--make-bed \
--out $outdir/target_variants


# calculate LD:
plink \
--bfile $outdir/target_variants \
--r2 square \
--out $outdir/target_variants


