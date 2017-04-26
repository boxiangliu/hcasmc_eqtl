plink --vcf ../../shared/1000genomes/phase3v5a/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep-allele-order --r2 --ld rs3131972 rs3131971

# --ld rs3131972 rs3131971:

#    R-sq = 0.617604       D' = 1

#    Haplotype     Frequency    Expectation under LE
#    ---------     ---------    --------------------
#           GC      0.653355                0.492103
#           AC      0.099840                0.261092
#           GT      0                       0.161251
#           AT      0.246805                0.085554

#    In phase alleles are GC/AT

# --r2 to plink.ld ... done.    
bcftools view ../../shared/1000genomes/phase3v5a/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 1:1-1752721 > ../processed_data/gwas_atacseq_overlap/tmp/test.vcf
bcftools view ../../shared/1000genomes/phase3v5a/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 1:1-1752721 > ../processed_data/gwas_atacseq_overlap/tmp/test.vcf

../../tools/bcftools/bcftools view ../../shared/1000genomes/phase3v5a/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 22:28698027 | ../../tools/bcftools/bcftools annotate  --set-id +'%CHROM:%POS' | ../../tools/bcftools/bcftools norm -d 'any'

plink --vcf ../processed_data/gwas_atacseq_overlap/tmp/test.vcf --keep-allele-order --r2 --ld-snps rs3131972,rs3131971 --ld-window-kb 1000 --ld-window-r2 0.5 --out ../processed_data/gwas_atacseq_overlap/tmp/test
vcftools --gzvcf ../processed_data/gwas_atacseq_overlap/tmp/chr22.vcf.gz --geno-r2-positions ../processed_data/gwas_atacseq_overlap/tmp/positions.chr22.tsv --ld-window-bp 1000000 --min-r2 0.7 --out ../processed_data/gwas_atacseq_overlap/tmp/ld0.7.vcftools