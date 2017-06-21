out_dir=../processed_data/finemap/finemap/preprocess
tmp_dir=../processed_data/finemap/finemap/tmp
[[ ! -d $tmp_dir ]] && mkdir -p $tmp_dir
[[ ! -d $out_dir ]] && mkdir -p $out_dir


# Subset to EUR and Howson variants: 
eur_sample_fn=$tmp_dir/eur_sample.txt
region_fn=$tmp_dir/howson.chrpos

grep EUR ../../shared/1000genomes/phase3v5a/integrated_call_samples_v3.20130502.ALL.panel | cut -f1 > $eur_sample_fn
tail -n +2 ../data/gwas/howson_2017/Howson-JMM_CHD_Mixed_2017.txt | cut -f2 | sed "s/chr//" | sed -r "s/:/\t/g" > $region_fn

for i in `seq 22`; do
../../tools/bcftools/bcftools view -S $eur_sample_fn \
-R $tmp_dir/howson.chrpos \
-v snps,indels,mnps \
../../shared/1000genomes/phase3v5a/ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
-Ou | ../../tools/bcftools/bcftools annotate --set-id '%CHROM:%POS' -Oz > $tmp_dir/chr$i.vcf.gz &
done

bcftools concat $tmp_dir/chr{1..22}.vcf.gz | bgzip > $tmp_dir/all.vcf.gz
rm $tmp_dir/chr{1..22}.vcf.gz


# Normalize and left align: 
bcftools norm --multiallelics -both --check-ref ws -f /srv/persistent/bliu2/shared/genomes/GRCh37/hs37d5.fa $tmp_dir/all.vcf.gz -Oz > $tmp_dir/all.norm.vcf.gz


echo -e "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > $tmp_dir/Howson.vcf
tail -n +2 ../data/gwas/howson_2017/Howson-JMM_CHD_Mixed_2017.txt | \
sed "s/chr//g" | awk 'BEGIN{FS="\t";OFS="\t"}{split($2,a,":"); print a[1],a[2],$1,$3,$4,".",".","."}' >> $tmp_dir/Howson.vcf
bcftools norm --multiallelics -both --check-ref ws -f /srv/persistent/bliu2/shared/genomes/GRCh37/hs37d5.fa $tmp_dir/Howson.vcf > $tmp_dir/Howson.norm.vcf


# Fix strand flips, reorder ref/alt, remove duplicate lines:
Rscript finemap/finemap/normalize_ref_alt.R


region_fn=$tmp_dir/howson.snpid
tail -n +2 ../data/gwas/howson_2017/Howson-JMM_CHD_Mixed_2017.norm.txt | cut -f13 | sed "s/chr//" > $region_fn
echo '#' >> $region_fn
../../tools/bcftools/bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz $tmp_dir/all.norm.vcf.gz > $tmp_dir/all.norm.id.vcf.gz


plink --extract $region_fn --vcf $tmp_dir/all.norm.id.vcf.gz --keep-allele-order --recode vcf --out $tmp_dir/all.final
rm $tmp_dir/all.final.{log,nosex}

bcftools view -H $tmp_dir/all.final.vcf | awk '{print "chr"$3}' > $tmp_dir/howson_in_1000g.snpid
grep -e rsid -wFf $tmp_dir/howson_in_1000g.snpid ../data/gwas/howson_2017/Howson-JMM_CHD_Mixed_2017.norm.txt > ../data/gwas/howson_2017/Howson-JMM_CHD_Mixed_2017.norm.in1kgp3.txt

ln $tmp_dir/all.final.vcf $out_dir

rm -r $tmp_dir
