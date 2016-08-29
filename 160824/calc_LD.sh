# algorithm: 
# subset to 100-200 top hits.

# for each gene:
#	determine chromosome
#	subset to relavant variants
# 	calculate LD
#	run eCavier

# filter: 
# for each gene: 
# 	read gwas zscore
#	determine whether at least 1 zscore cross the threshold -> bool A
# 	read eqtl zscore 
#	determine whether at least 1 zscore cross the threshold -> bool B
# 	gene is protein coding or lincRNA -> bool C
#	output gene name if A&B&C



# determin chromsome: 
# for each gene: 
#	read the first line of the zscore file
#	parse the first line, find the chromsome
# 	output <gene_name> <chromosome> 



# check if all snps have variant id
# subset to relevant variants: 
cat $snp_list | awk 'BEGIN{FS="_";OFS="\t"}{print $1,$2}' > ${snp_list}.tmp
cat $eur_sample | awk '{print $1}' > $eur_sample.tmp
bcftools view -R ${snp_list}.tmp -S $eur_sample.tmp -Ov -o /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160824/ld/ENSG00000000419.8.vcf $in_file


# calculate LD
in_file=/srv/persistent/bliu2/shared/1000genomes/phase1v3/ALL.chr20.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.rsid.vcf.gz
snp_list=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160824/eCAVIAR_input/ENSG00000000419.8.id
eur_sample=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160824/eur_samples.txt
plink --vcf /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160824/ld/ENSG00000000419.8.vcf \
	--make-bed \
	--double-id \
	--extract $snp_list \



plink -bfile plink --r square

