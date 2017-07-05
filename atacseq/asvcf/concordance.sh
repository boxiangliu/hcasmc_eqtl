out_dir=../processed_data/atacseq/asvcf/concordance/
mkdir -p $out_dir

bcftools view -m2 -M2 -Ou ../data/joint3/asvcf_sid_atac/phased_and_imputed.chr1.rename.dr2.hwe.indellt51.atacsample.hg19.vcf.new.gz | bcftools query -f '%CHROM\_%POS\_%REF\_%ALT[\t%DS]\n'  > $out_dir/genotype.txt
bcftools view -m2 -M2 -Ou ../data/joint3/asvcf_sid_atac/phased_and_imputed.chr1.rename.dr2.hwe.indellt51.atacsample.hg19.vcf.new.gz | bcftools query -f '%CHROM\_%POS\_%REF\_%ALT[\t%AS]\n'  > $out_dir/ase.txt

Rscript atacseq/asvcf/concordance.R