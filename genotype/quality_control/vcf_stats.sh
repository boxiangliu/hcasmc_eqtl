# Bcftools stats:
bcftools stats ../data/joint2/recalibrated_variants.pass.vcf.gz > ../processed_data/160514_Ts_Tv_ratio/stats.txt

# calculate Ti/Tv ratio:
plot-vcfstats --no-PDF -p ../figures/genotype/joint_genotyping/quality_control/vcf_stats ../processed_data/genotype/quality_control/stats.txt
