in_dir=/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/asvcf_sid/
atac_dir=/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/asvcf_sid_atac/

# Subset to ATACseq sample:
parallel -j 22 \
bcftools view -s 1346,1508,1522,200212,2305,2356,2510,2989 \
-Ou $in_dir/phased_and_imputed.chr{}.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz "|" bcftools annotate -x FORMAT/AS -Oz \
-o $atac_dir/phased_and_imputed.chr{}.rename.dr2.hwe.indellt51.atacsample.hg19.vcf.gz ::: {1..22}

# Index:
parallel -j 22 \
tabix -p vcf $atac_dir/phased_and_imputed.chr{}.rename.dr2.hwe.indellt51.atacsample.hg19.vcf.gz ::: {1..22}


# Count AS reads:
parallel -j 22 \
bash /srv/persistent/bliu2/tools/rasqual/src/ASVCF/createASVCF.sh \
/srv/persistent/bliu2/HCASMC_eQTL/scripts/atacseq/asvcf/bam.list.txt \
$atac_dir/phased_and_imputed.chr{}.rename.dr2.hwe.indellt51.atacsample.hg19.vcf.gz ::: {1..22}


parallel -j 22 \
tabix -p vcf $atac_dir/phased_and_imputed.chr{}.rename.dr2.hwe.indellt51.atacsample.hg19.vcf.new.gz ::: {1..22}
