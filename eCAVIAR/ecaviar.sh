# command args: 
in_dir=$1
out_dir=$2
echo "input directory: $in_dir"
echo "output directory: $out_dir"

# in_dir=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/eCAVIAR/eCAVIAR_input/
# out_dir=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/eCAVIAR/eCAVIAR_output/


# path to 1000 Genome European sample IDs:
eur_sample=/srv/persistent/bliu2/HCASMC_eQTL/processed_data//eCAVIAR/eur_samples.txt


# loop over every SNP list: 
for snp_list in $(ls $in_dir/*.id); do
	base=$(basename $snp_list)

	# if [[ ! "$base" == "ENSG00000031698.8.id" ]]; then continue; fi
	echo "SNP list: $snp_list" 
	i=$(head -n1 $snp_list | awk 'BEGIN{FS="_"}{print $1}')
	in_file=/srv/persistent/bliu2/shared/1000genomes/phase1v3/ALL.chr$i.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.rsid.vcf.gz


	# split snp into chr and pos: 
	cat $snp_list | awk 'BEGIN{FS="_";OFS="\t"}{print $1,$2}' > ${snp_list}.tmp
	cat $eur_sample | awk '{print $1}' > $eur_sample.tmp


	# subset to relevant variants: 
	echo "subsetting vcf..."
	bcftools view -R ${snp_list}.tmp -S $eur_sample.tmp -Ov -o ${snp_list/id/vcf} $in_file


	# generate bed bim fam:
	echo "generating bed, bim and fam files..."
	plink --vcf ${snp_list/id/vcf} \
		--keep-allele-order \
		--make-bed \
		--double-id \
		--extract $snp_list \
		--out ${snp_list/.id/}


	# calculate LD:
	echo "calculating LD..."
	plink -bfile ${snp_list/.id/} \
		--r square \
		--out ${snp_list/.id/}


	# run eCAVIAR:
	echo "running eCAVIAR..."
	base=$(basename $snp_list)
	if [[ ! -d $out_dir ]]; then mkdir -p $out_dir; fi
	/srv/persistent/bliu2/tools/caviar/CAVIAR-C++/eCAVIAR \
		-o $out_dir/${base/.id/.ecaviar} \
		-l ${snp_list/id/ld} \
		-l ${snp_list/id/ld} \
		-z ${snp_list/id/gwas.zscore} \
		-z ${snp_list/id/eqtl.zscore} \
		-r 0.95
done



