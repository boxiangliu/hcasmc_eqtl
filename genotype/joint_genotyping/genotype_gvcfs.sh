in_dir=../data/joint2 

# Genotype GVCF for each chromosome: 
for chr in $(seq 1 22) X Y M; do 
	bash genotype/joint_genotyping/genotype_gvcfs.core.sh chr$chr 2> ../logs/genotype/joint_genotyping/genotype_gvcfs.$chr.$dat.log &
done

mkdir by_chrom
mv $in_dir/raw_variants.chr*.vcf* $in_dir/by_chrom

# Concatenate all chromosomes: 
bcftools concat -Oz -o $in_dir/raw_variants.vcf.gz $in_dir/by_chrom/raw_variants.chr{M,{1..22},X,Y}.vcf.gz


